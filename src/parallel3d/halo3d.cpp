// Function bodies for parallel3d/halo3d.hpp.

#include "parallel3d/halo3d.hpp"

#include <cstddef>

namespace {

// One tag per direction of travel, as in the 2D exchange.
enum Tag { moving_north = 11, moving_south = 12, moving_west = 13, moving_east = 14 };

static_assert(halo == 1, "Halo3D exchanges one layer of ghost cells");

}  // namespace

Halo3D::Halo3D(const ProcessLayout& layout, int nz, int ny, int nx)
    : comm_(layout.comm),
      north_(layout.north),
      south_(layout.south),
      west_(layout.west),
      east_(layout.east),
      nz_(nz),
      ny_(ny),
      nx_(nx),
      rows_(MPI_DATATYPE_NULL),
      send_west_(static_cast<std::size_t>(nz) * ny),
      send_east_(static_cast<std::size_t>(nz) * ny),
      recv_west_(static_cast<std::size_t>(nz) * ny),
      recv_east_(static_cast<std::size_t>(nz) * ny) {
    const std::size_t layer = static_cast<std::size_t>(ny_ + 2 * halo) * static_cast<std::size_t>(nx_ + 2 * halo);
    mpi_check(MPI_Type_vector(nz_, nx_, static_cast<int>(layer), MPI_DOUBLE, &rows_), "MPI_Type_vector (3D rows)");
    mpi_check(MPI_Type_commit(&rows_), "MPI_Type_commit (3D rows)");
}

Halo3D::~Halo3D() {
    // No mpi_check: a destructor must not throw.
    if (rows_ != MPI_DATATYPE_NULL) MPI_Type_free(&rows_);
}

// The west buffer holds my column 1 of every row of every layer, which
// becomes my west neighbour's ghost column nx+1; the east buffer holds my
// column nx, which becomes my east neighbour's ghost column 0.
void Halo3D::pack_columns(const Field3D& c) {
    for (int k = 1; k <= nz_; ++k) {
        for (int i = 1; i <= ny_; ++i) {
            const std::size_t at = static_cast<std::size_t>(k - 1) * ny_ + (i - 1);
            const double* row = c.row(k, i);
            send_west_[at] = row[1];
            send_east_[at] = row[nx_];
        }
    }
}

void Halo3D::unpack_columns(Field3D& c) {
    for (int k = 1; k <= nz_; ++k) {
        for (int i = 1; i <= ny_; ++i) {
            const std::size_t at = static_cast<std::size_t>(k - 1) * ny_ + (i - 1);
            double* row = c.row(k, i);
            if (west_ != MPI_PROC_NULL) row[0] = recv_west_[at];
            if (east_ != MPI_PROC_NULL) row[nx_ + 1] = recv_east_[at];
        }
    }
}

void Halo3D::exchange(Field3D& c) {
    const int count = nz_ * ny_;
    // My top rows go north; my bottom ghost rows come from the south. Then the reverse.
    mpi_check(MPI_Sendrecv(&c(1, 1, 1), 1, rows_, north_, moving_north, &c(1, ny_ + 1, 1), 1, rows_, south_,
                           moving_north, comm_, MPI_STATUS_IGNORE),
              "3D halo north/south");
    mpi_check(MPI_Sendrecv(&c(1, ny_, 1), 1, rows_, south_, moving_south, &c(1, 0, 1), 1, rows_, north_, moving_south,
                           comm_, MPI_STATUS_IGNORE),
              "3D halo south/north");
    pack_columns(c);
    mpi_check(MPI_Sendrecv(send_west_.data(), count, MPI_DOUBLE, west_, moving_west, recv_east_.data(), count,
                           MPI_DOUBLE, east_, moving_west, comm_, MPI_STATUS_IGNORE),
              "3D halo west/east");
    mpi_check(MPI_Sendrecv(send_east_.data(), count, MPI_DOUBLE, east_, moving_east, recv_west_.data(), count,
                           MPI_DOUBLE, west_, moving_east, comm_, MPI_STATUS_IGNORE),
              "3D halo east/west");
    unpack_columns(c);
}

void Halo3D::start(Field3D& c) {
    const int count = nz_ * ny_;
    pending_ = 0;
    mpi_check(MPI_Irecv(&c(1, 0, 1), 1, rows_, north_, moving_south, comm_, &requests_[pending_++]), "3D Irecv north");
    mpi_check(MPI_Irecv(&c(1, ny_ + 1, 1), 1, rows_, south_, moving_north, comm_, &requests_[pending_++]),
              "3D Irecv south");
    mpi_check(MPI_Irecv(recv_west_.data(), count, MPI_DOUBLE, west_, moving_east, comm_, &requests_[pending_++]),
              "3D Irecv west");
    mpi_check(MPI_Irecv(recv_east_.data(), count, MPI_DOUBLE, east_, moving_west, comm_, &requests_[pending_++]),
              "3D Irecv east");
    pack_columns(c);
    mpi_check(MPI_Isend(&c(1, 1, 1), 1, rows_, north_, moving_north, comm_, &requests_[pending_++]), "3D Isend north");
    mpi_check(MPI_Isend(&c(1, ny_, 1), 1, rows_, south_, moving_south, comm_, &requests_[pending_++]),
              "3D Isend south");
    mpi_check(MPI_Isend(send_west_.data(), count, MPI_DOUBLE, west_, moving_west, comm_, &requests_[pending_++]),
              "3D Isend west");
    mpi_check(MPI_Isend(send_east_.data(), count, MPI_DOUBLE, east_, moving_east, comm_, &requests_[pending_++]),
              "3D Isend east");
}

void Halo3D::finish(Field3D& c) {
    mpi_check(MPI_Waitall(pending_, requests_, MPI_STATUSES_IGNORE), "3D MPI_Waitall halo");
    pending_ = 0;
    unpack_columns(c);
}
