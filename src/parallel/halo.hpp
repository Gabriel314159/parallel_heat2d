#pragma once
#include "parallel/decomposition.hpp"
#include <array>

namespace heat2d {
class Halo {
public:
    explicit Halo(const Domain& domain);
    ~Halo();
    Halo(const Halo&)=delete;
    Halo& operator=(const Halo&)=delete;
    void blocking(Grid& g);
    void start(Grid& g);
    void finish();
private:
    const Domain& d_;
    MPI_Datatype column_=MPI_DATATYPE_NULL;
    std::array<MPI_Request,8> requests_{};
    bool active_=false;
};
}
