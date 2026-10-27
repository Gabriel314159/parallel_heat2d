// Function bodies for parallel/halo.hpp.

#include "parallel/halo.hpp"

#include <cstddef>

namespace {

// One tag per direction of travel, so the four messages a rank receives in a
// step can never be confused with each other.
enum Tag { moving_north = 1, moving_south = 2, moving_west = 3, moving_east = 4 };

}  // namespace

HaloExchange::HaloExchange(const ProcessLayout& layout, int ny, int nx)
    : comm_(layout.comm),
      north_(layout.north),
      south_(layout.south),
      west_(layout.west),
      east_(layout.east),
      ny_(ny),
      nx_(nx),
      rows_(MPI_DATATYPE_NULL),
      send_west_(static_cast<std::size_t>(halo) * ny),
      send_east_(static_cast<std::size_t>(halo) * ny),
      recv_west_(static_cast<std::size_t>(halo) * ny),
      recv_east_(static_cast<std::size_t>(halo) * ny) {
    MPI_Type_vector(halo, nx_, nx_ + 2 * halo, MPI_DOUBLE, &rows_);
    MPI_Type_commit(&rows_);
}

HaloExchange::~HaloExchange() {
    if (rows_ != MPI_DATATYPE_NULL) MPI_Type_free(&rows_);
}

// Layer k = 1..halo counts inwards from the edge. The west buffer holds my
// columns 1..halo, which become my west neighbour's ghost columns nx+1..nx+halo;
// the east buffer holds columns nx, nx-1, ..., which become my east
// neighbour's ghost columns 0, -1, ....
void HaloExchange::pack_columns(const Field& c) {
    for (int k = 1; k <= halo; ++k) {
        for (int i = 1; i <= ny_; ++i) {
            const std::size_t at = static_cast<std::size_t>(k - 1) * ny_ + (i - 1);
            send_west_[at] = c(i, k);
            send_east_[at] = c(i, nx_ + 1 - k);
        }
    }
}

void HaloExchange::unpack_columns(Field& c) {
    for (int k = 1; k <= halo; ++k) {
        for (int i = 1; i <= ny_; ++i) {
            const std::size_t at = static_cast<std::size_t>(k - 1) * ny_ + (i - 1);
            if (west_ != MPI_PROC_NULL) c(i, 1 - k) = recv_west_[at];
            if (east_ != MPI_PROC_NULL) c(i, nx_ + k) = recv_east_[at];
        }
    }
}

void HaloExchange::exchange(Field& c) {
    const int count = halo * ny_;
    // My top rows go north; my bottom ghost rows come from the south. Then the reverse.
    MPI_Sendrecv(&c(1, 1), 1, rows_, north_, moving_north, &c(ny_ + 1, 1), 1, rows_, south_, moving_north, comm_,
                 MPI_STATUS_IGNORE);
    MPI_Sendrecv(&c(ny_ - halo + 1, 1), 1, rows_, south_, moving_south, &c(1 - halo, 1), 1, rows_, north_,
                 moving_south, comm_, MPI_STATUS_IGNORE);
    pack_columns(c);
    MPI_Sendrecv(send_west_.data(), count, MPI_DOUBLE, west_, moving_west, recv_east_.data(), count, MPI_DOUBLE,
                 east_, moving_west, comm_, MPI_STATUS_IGNORE);
    MPI_Sendrecv(send_east_.data(), count, MPI_DOUBLE, east_, moving_east, recv_west_.data(), count, MPI_DOUBLE,
                 west_, moving_east, comm_, MPI_STATUS_IGNORE);
    unpack_columns(c);
}

void HaloExchange::start(Field& c) {
    const int count = halo * ny_;
    pending_ = 0;
    MPI_Irecv(&c(1 - halo, 1), 1, rows_, north_, moving_south, comm_, &requests_[pending_++]);
    MPI_Irecv(&c(ny_ + 1, 1), 1, rows_, south_, moving_north, comm_, &requests_[pending_++]);
    MPI_Irecv(recv_west_.data(), count, MPI_DOUBLE, west_, moving_east, comm_, &requests_[pending_++]);
    MPI_Irecv(recv_east_.data(), count, MPI_DOUBLE, east_, moving_west, comm_, &requests_[pending_++]);
    pack_columns(c);
    MPI_Isend(&c(1, 1), 1, rows_, north_, moving_north, comm_, &requests_[pending_++]);
    MPI_Isend(&c(ny_ - halo + 1, 1), 1, rows_, south_, moving_south, comm_, &requests_[pending_++]);
    MPI_Isend(send_west_.data(), count, MPI_DOUBLE, west_, moving_west, comm_, &requests_[pending_++]);
    MPI_Isend(send_east_.data(), count, MPI_DOUBLE, east_, moving_east, comm_, &requests_[pending_++]);
}

void HaloExchange::finish(Field& c) {
    MPI_Waitall(pending_, requests_, MPI_STATUSES_IGNORE);
    pending_ = 0;
    unpack_columns(c);
}
