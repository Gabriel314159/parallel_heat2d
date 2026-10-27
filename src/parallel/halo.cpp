#include "parallel/halo.hpp"

namespace heat2d {
Halo::Halo(const Domain& domain): d_(domain) {
    mpi_check(MPI_Type_vector(d_.ny,1,d_.nx+2,MPI_DOUBLE,&column_),"MPI_Type_vector");
    mpi_check(MPI_Type_commit(&column_),"MPI_Type_commit");
}
Halo::~Halo() {
    if (active_) MPI_Waitall(8,requests_.data(),MPI_STATUSES_IGNORE);
    if (column_!=MPI_DATATYPE_NULL) MPI_Type_free(&column_);
}
void Halo::blocking(Grid& g) {
    mpi_check(MPI_Sendrecv(g.row(1)+1,g.nx,MPI_DOUBLE,d_.north,100,
        g.row(g.ny+1)+1,g.nx,MPI_DOUBLE,d_.south,100,d_.comm,MPI_STATUS_IGNORE),"halo north/south");
    mpi_check(MPI_Sendrecv(g.row(g.ny)+1,g.nx,MPI_DOUBLE,d_.south,101,
        g.row(0)+1,g.nx,MPI_DOUBLE,d_.north,101,d_.comm,MPI_STATUS_IGNORE),"halo south/north");
    mpi_check(MPI_Sendrecv(&g(1,1),1,column_,d_.west,102,
        &g(1,g.nx+1),1,column_,d_.east,102,d_.comm,MPI_STATUS_IGNORE),"halo west/east");
    mpi_check(MPI_Sendrecv(&g(1,g.nx),1,column_,d_.east,103,
        &g(1,0),1,column_,d_.west,103,d_.comm,MPI_STATUS_IGNORE),"halo east/west");
}
void Halo::start(Grid& g) {
    if (active_) throw std::logic_error("Halo exchange is already active.");
    requests_.fill(MPI_REQUEST_NULL);
    active_=true;
    mpi_check(MPI_Irecv(g.row(0)+1,g.nx,MPI_DOUBLE,d_.north,101,d_.comm,&requests_[0]),"Irecv north");
    mpi_check(MPI_Irecv(g.row(g.ny+1)+1,g.nx,MPI_DOUBLE,d_.south,100,d_.comm,&requests_[1]),"Irecv south");
    mpi_check(MPI_Irecv(&g(1,0),1,column_,d_.west,103,d_.comm,&requests_[2]),"Irecv west");
    mpi_check(MPI_Irecv(&g(1,g.nx+1),1,column_,d_.east,102,d_.comm,&requests_[3]),"Irecv east");
    mpi_check(MPI_Isend(g.row(1)+1,g.nx,MPI_DOUBLE,d_.north,100,d_.comm,&requests_[4]),"Isend north");
    mpi_check(MPI_Isend(g.row(g.ny)+1,g.nx,MPI_DOUBLE,d_.south,101,d_.comm,&requests_[5]),"Isend south");
    mpi_check(MPI_Isend(&g(1,1),1,column_,d_.west,102,d_.comm,&requests_[6]),"Isend west");
    mpi_check(MPI_Isend(&g(1,g.nx),1,column_,d_.east,103,d_.comm,&requests_[7]),"Isend east");
}
void Halo::finish() {
    if (active_) {
        mpi_check(MPI_Waitall(8,requests_.data(),MPI_STATUSES_IGNORE),"MPI_Waitall halo");
        active_=false;
    }
}
}
