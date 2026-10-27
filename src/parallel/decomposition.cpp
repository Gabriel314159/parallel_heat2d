#include "parallel/decomposition.hpp"
#include "common/physics.hpp"
#include <algorithm>
#include <limits>
#include <vector>

namespace heat2d {
namespace {
int extent(int n, int blocks, int coord) { return n/blocks + (coord<n%blocks ? 1:0); }
int offset(int n, int blocks, int coord) { return coord*(n/blocks)+std::min(coord,n%blocks); }
}
Domain::Domain(const Config& c, MPI_Comm world) {
    mpi_check(MPI_Comm_size(world,&size),"MPI_Comm_size");
    if (c.decomposition=="rows") {
        py=size; px=1;
        if ((c.py && c.py!=py) || (c.px && c.px!=px))
            throw std::invalid_argument("Row decomposition requires py=ranks and px=1.");
    } else if (c.py || c.px) {
        py=c.py; px=c.px;
        if (!py) { if (size%px) throw std::invalid_argument("px must divide ranks."); py=size/px; }
        if (!px) { if (size%py) throw std::invalid_argument("py must divide ranks."); px=size/py; }
    } else {
        // Minimize approximate per-rank perimeter, respecting thin global grids.
        double best=std::numeric_limits<double>::infinity();
        py=0; px=0;
        for (int y=1; y<=size; ++y) {
            if (size%y) continue;
            const int x=size/y;
            if (y>c.ny || x>c.nx) continue;
            const double cost=double(c.ny)/y+double(c.nx)/x;
            if (cost<best) { best=cost; py=y; px=x; }
        }
    }
    if (py<1 || px<1 || static_cast<long long>(py)*px!=size || py>c.ny || px>c.nx)
        throw std::invalid_argument("Process grid must cover all ranks and give every rank a nonempty block.");
    int dims[]={py,px}, periodic[]={0,0};
    mpi_check(MPI_Cart_create(world,2,dims,periodic,0,&comm),"MPI_Cart_create");
    mpi_check(MPI_Comm_set_errhandler(comm,MPI_ERRORS_RETURN),"MPI_Comm_set_errhandler");
    mpi_check(MPI_Comm_rank(comm,&rank),"MPI_Comm_rank");
    int coords[2]; mpi_check(MPI_Cart_coords(comm,rank,2,coords),"MPI_Cart_coords");
    ny=extent(c.ny,py,coords[0]); nx=extent(c.nx,px,coords[1]);
    oy=offset(c.ny,py,coords[0]); ox=offset(c.nx,px,coords[1]);
    mpi_check(MPI_Cart_shift(comm,0,1,&north,&south),"MPI_Cart_shift rows");
    mpi_check(MPI_Cart_shift(comm,1,1,&west,&east),"MPI_Cart_shift columns");
}
Domain::~Domain() { if (comm!=MPI_COMM_NULL) MPI_Comm_free(&comm); }

std::unique_ptr<Grid> gather_global(const Grid& local, const Config& c, const Domain& d) {
    const long long global_count=static_cast<long long>(c.nx)*c.ny;
    if (global_count>std::numeric_limits<int>::max())
        throw std::runtime_error("CSV gathering uses MPI int counts; use a smaller output grid.");
    const int count=local.nx*local.ny;
    std::vector<double> packed(static_cast<std::size_t>(count));
    for (int i=1; i<=local.ny; ++i)
        std::copy(local.row(i)+1,local.row(i)+1+local.nx,packed.begin()+(i-1)*local.nx);
    int meta[]={d.oy,d.ox,d.ny,d.nx};
    std::vector<int> metadata(d.rank==0 ? 4*d.size:0);
    mpi_check(MPI_Gather(meta,4,MPI_INT,metadata.data(),4,MPI_INT,0,d.comm),"MPI_Gather metadata");
    std::vector<int> counts, displacements;
    std::vector<double> values;
    if (d.rank==0) {
        counts.resize(d.size); displacements.resize(d.size);
        int position=0;
        for (int r=0; r<d.size; ++r) {
            counts[r]=metadata[4*r+2]*metadata[4*r+3];
            displacements[r]=position; position+=counts[r];
        }
        values.resize(static_cast<std::size_t>(position));
    }
    mpi_check(MPI_Gatherv(packed.data(),count,MPI_DOUBLE,values.data(),counts.data(),
                         displacements.data(),MPI_DOUBLE,0,d.comm),"MPI_Gatherv field");
    if (d.rank!=0) return nullptr;
    auto global=std::make_unique<Grid>(c.nx,c.ny);
    initialize(*global,c);
    for (int r=0; r<d.size; ++r) {
        const int oy=metadata[4*r], ox=metadata[4*r+1];
        const int ny=metadata[4*r+2], nx=metadata[4*r+3];
        for (int i=0; i<ny; ++i)
            std::copy_n(values.data()+displacements[r]+i*nx,nx,global->row(oy+i+1)+ox+1);
    }
    return global;
}
}
