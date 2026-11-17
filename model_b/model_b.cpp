// Model B: staggered 3D incompressible flow, pressure projection and passive pigment.
// All coordinate axes and array indices increase together. MPI owns each high face
// exactly once, at the cell on its low side; the low boundary faces are zero ghosts.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
using Clock = std::chrono::steady_clock;
constexpr double pi = 3.14159265358979323846;
double elapsed(Clock::time_point t) { return std::chrono::duration<double>(Clock::now()-t).count(); }

struct Options {
    std::array<int,3> n{32,32,16}, dims{0,0,0};
    std::array<double,3> length{1,1,0.5};
    int steps=100, spinup=0, max_iter=10000, progress=0;
    double nu=1e-3, kappa=1e-4, force=0.8, tol=1e-10, atol=1e-13;
    double dt=0, speed_bound=2, wind_alpha=0, wind_omega=1;
    std::string mode="simulate", profile="surface", initial="blob", decomp="blocks";
    std::string output, trace;
};
Options parse(int argc,char** argv) {
    Options o;
    for(int a=1;a<argc;++a) {
        std::string k=argv[a];
        if(k=="--help") {
            std::cout << "--nx N --ny N --nz N --lx L --ly L --lz L --steps N --spinup N\n"
                      << "--nu V --kappa K --force F --dt DT --speed-bound U --tol R --atol A\n"
                      << "--profile surface|uniform --initial blob|constant|cosine\n"
                      << "--decomp slabs|pencils|blocks --px P --py P --pz P\n"
                      << "--mode simulate|poisson-discrete|poisson-continuous|projection|momentum-test\n"
                      << "--output final.csv --trace trace.csv --max-iter N --progress N\n";
            o.steps=-1; return o;
        }
        if(a+1==argc) throw std::runtime_error("missing value for "+k);
        std::string v=argv[++a];
        if(k=="--nx") o.n[0]=std::stoi(v);
        else if(k=="--ny") o.n[1]=std::stoi(v);
        else if(k=="--nz") o.n[2]=std::stoi(v);
        else if(k=="--lx") o.length[0]=std::stod(v);
        else if(k=="--ly") o.length[1]=std::stod(v);
        else if(k=="--lz") o.length[2]=std::stod(v);
        else if(k=="--px") o.dims[0]=std::stoi(v);
        else if(k=="--py") o.dims[1]=std::stoi(v);
        else if(k=="--pz") o.dims[2]=std::stoi(v);
        else if(k=="--steps") o.steps=std::stoi(v);
        else if(k=="--spinup") o.spinup=std::stoi(v);
        else if(k=="--max-iter") o.max_iter=std::stoi(v);
        else if(k=="--progress") o.progress=std::stoi(v);
        else if(k=="--nu") o.nu=std::stod(v);
        else if(k=="--kappa") o.kappa=std::stod(v);
        else if(k=="--force") o.force=std::stod(v);
        else if(k=="--tol") o.tol=std::stod(v);
        else if(k=="--atol") o.atol=std::stod(v);
        else if(k=="--dt") o.dt=std::stod(v);
        else if(k=="--speed-bound") o.speed_bound=std::stod(v);
        else if(k=="--wind-alpha") o.wind_alpha=std::stod(v);
        else if(k=="--wind-omega") o.wind_omega=std::stod(v);
        else if(k=="--mode") o.mode=v;
        else if(k=="--profile") o.profile=v;
        else if(k=="--initial") o.initial=v;
        else if(k=="--decomp") o.decomp=v;
        else if(k=="--output") o.output=v;
        else if(k=="--trace") o.trace=v;
        else throw std::runtime_error("unknown option "+k);
    }
    for(int d=0;d<3;++d)
        if(o.n[d]<2 || o.dims[d]<0 || !std::isfinite(o.length[d]) || o.length[d]<=0)
            throw std::runtime_error("invalid grid / length / process grid");
    for(double x:{o.nu,o.kappa,o.force,o.dt,o.speed_bound,o.tol,o.atol,o.wind_alpha,o.wind_omega})
        if(!std::isfinite(x)) throw std::runtime_error("non-finite parameter");
    if(o.steps<1 || o.spinup<0 || o.max_iter<1 || o.progress<0 || o.nu<0 || o.kappa<0 || o.dt<0 ||
       o.speed_bound<=0 || o.tol<=0 || o.tol>=1 || o.atol<=0 ||
       (o.mode!="simulate" && o.mode!="poisson-discrete" && o.mode!="poisson-continuous" && o.mode!="projection" && o.mode!="momentum-test") ||
       (o.profile!="surface" && o.profile!="uniform") ||
       (o.initial!="blob" && o.initial!="constant" && o.initial!="cosine") ||
       (o.decomp!="slabs" && o.decomp!="pencils" && o.decomp!="blocks"))
        throw std::runtime_error("invalid physical or run parameter");
    return o;
}
struct Domain {
    int rank=0,ranks=1,nodes=1;
    std::array<int,3> n{},start{},dims{1,1,1},coord{},global{};
    std::array<double,3> h{};
    double count=0,comm_s=0,reduce_s=0;
#ifdef USE_MPI
    MPI_Comm cart=MPI_COMM_NULL;
    int peer[3][2]{};
#endif
    explicit Domain(const Options& o):global(o.n) {
#ifdef USE_MPI
        MPI_Comm_size(MPI_COMM_WORLD,&ranks);
#endif
        // Search feasible factorizations, including uneven and very thin blocks.
        double best=std::numeric_limits<double>::infinity();
        for(int px=1;px<=ranks;++px) if(ranks%px==0)
        for(int py=1;py<=ranks/px;++py) if((ranks/px)%py==0) {
            int pz=ranks/px/py;
            std::array<int,3> p{px,py,pz};
            bool valid=true;
            for(int d=0;d<3;++d) valid=valid && p[d]<=o.n[d] && (!o.dims[d] || o.dims[d]==p[d]);
            if(o.decomp=="slabs") valid=valid && px==1 && pz==1;
            if(o.decomp=="pencils") valid=valid && pz==1;
            if(!valid) continue;
            double cost=0;
            for(int d=0;d<3;++d) cost+=(p[d]-1.0)/o.n[d];
            if(cost<best) { best=cost; dims=p; }
        }
        if(!std::isfinite(best)) throw std::runtime_error("no feasible MPI process grid");
#ifdef USE_MPI
        int periods[3]={0,0,0};
        MPI_Cart_create(MPI_COMM_WORLD,3,dims.data(),periods,0,&cart);
        MPI_Comm_rank(cart,&rank);
        MPI_Cart_coords(cart,rank,3,coord.data());
        for(int d=0;d<3;++d) MPI_Cart_shift(cart,d,1,&peer[d][0],&peer[d][1]);
        MPI_Comm local;
        MPI_Comm_split_type(cart,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&local);
        int lr=0;
        MPI_Comm_rank(local,&lr);
        nodes=static_cast<int>(sum(lr==0 ? 1.0:0.0));
        MPI_Comm_free(&local);
#endif
        count=1.0*o.n[0]*o.n[1]*o.n[2];
        for(int d=0;d<3;++d) {
            n[d]=o.n[d]/dims[d]+(coord[d]<o.n[d]%dims[d]);
            start[d]=coord[d]*(o.n[d]/dims[d])+std::min(coord[d],o.n[d]%dims[d]);
            h[d]=o.length[d]/o.n[d];
        }
    }
    bool low(int d) const { return coord[d]==0; }
    bool high(int d) const { return coord[d]+1==dims[d]; }
    double sum(double v) {
#ifdef USE_MPI
        if(ranks>1) { auto t=Clock::now(); double z; MPI_Allreduce(&v,&z,1,MPI_DOUBLE,MPI_SUM,cart); reduce_s+=elapsed(t); return z; }
#endif
        return v;
    }
    double maximum(double v) {
#ifdef USE_MPI
        if(ranks>1) { auto t=Clock::now(); double z; MPI_Allreduce(&v,&z,1,MPI_DOUBLE,MPI_MAX,cart); reduce_s+=elapsed(t); return z; }
#endif
        return v;
    }
    void barrier() {
#ifdef USE_MPI
        MPI_Barrier(cart);
#endif
    }
    ~Domain() {
#ifdef USE_MPI
        if(cart!=MPI_COMM_NULL) MPI_Comm_free(&cart);
#endif
    }
};
struct Field {
    int nx,ny,nz;
    std::size_t sy,sz,size;
    std::unique_ptr<double[]> values;
    explicit Field(const Domain& d):nx(d.n[0]),ny(d.n[1]),nz(d.n[2]),sy(nx+2),sz(sy*(ny+2)),size(sz*(nz+2)),values(new double[size]) {
        #pragma omp parallel for schedule(static)
        for(std::size_t t=0;t<size;++t) values[t]=0;
    }
    std::size_t index(int i,int j,int k) const { return i+sy*j+sz*k; }
    double& operator()(int i,int j,int k) { return values[index(i,j,k)]; }
    double operator()(int i,int j,int k) const { return values[index(i,j,k)]; }
    double& at(const std::array<int,3>& q) { return (*this)(q[0],q[1],q[2]); }
    double at(const std::array<int,3>& q) const { return (*this)(q[0],q[1],q[2]); }
    void swap(Field& f) { values.swap(f.values); }
};
#define CELLS(f) for(int k=1;k<=(f).nz;++k) for(int j=1;j<=(f).ny;++j) for(int i=1;i<=(f).nx;++i)

// Sequential direction exchanges propagate edge/corner ghosts for momentum.
// Scalar Laplacian/transport only needs face ghosts, so it sends owned cross-sections.
void halo(Field& f,Domain& d,int component=-1,bool corners=false) {
    auto t=Clock::now();
    for(int a=0;a<3;++a) {
        int b=(a+1)%3,c=(a+2)%3;
        int lo=corners?0:1, nb=d.n[b]+(corners?2:0), nc=d.n[c]+(corners?2:0);
        std::size_t len=static_cast<std::size_t>(nb)*nc;
        if(len>static_cast<std::size_t>(std::numeric_limits<int>::max())) throw std::runtime_error("halo exceeds MPI int count");
        std::vector<double> sendlo(len),sendhi(len),recvlo(len),recvhi(len);
        std::size_t pos=0;
        for(int ic=lo;ic<lo+nc;++ic) for(int ib=lo;ib<lo+nb;++ib) {
            std::array<int,3> q{}; q[b]=ib; q[c]=ic;
            q[a]=1; sendlo[pos]=f.at(q);
            q[a]=d.n[a]; sendhi[pos++]=f.at(q);
        }
#ifdef USE_MPI
        MPI_Sendrecv(sendlo.data(),static_cast<int>(len),MPI_DOUBLE,d.peer[a][0],2*a,
                     recvhi.data(),static_cast<int>(len),MPI_DOUBLE,d.peer[a][1],2*a,d.cart,MPI_STATUS_IGNORE);
        MPI_Sendrecv(sendhi.data(),static_cast<int>(len),MPI_DOUBLE,d.peer[a][1],2*a+1,
                     recvlo.data(),static_cast<int>(len),MPI_DOUBLE,d.peer[a][0],2*a+1,d.cart,MPI_STATUS_IGNORE);
#endif
        pos=0;
        for(int ic=lo;ic<lo+nc;++ic) for(int ib=lo;ib<lo+nb;++ib) {
            std::array<int,3> q{}; q[b]=ib; q[c]=ic;
            q[a]=0;
            f.at(q)=d.low(a)?(component==a?0:sendlo[pos]):recvlo[pos];
            q[a]=d.n[a]+1;
            // Reflection of a normal component outside the high wall; the wall face itself is zero.
            if(d.high(a) && component==a) { auto r=q; r[a]=d.n[a]-1; f.at(q)=-f.at(r); }
            else f.at(q)=d.high(a)?sendhi[pos]:recvhi[pos];
            ++pos;
        }
    }
    d.comm_s+=elapsed(t);
}
double dot(const Field& a,const Field& b,Domain& d) {
    double s=0;
    #pragma omp parallel for collapse(3) reduction(+:s) schedule(static)
    CELLS(a) s+=a(i,j,k)*b(i,j,k);
    return d.sum(s);
}
double mean(Field& a,Domain& d) {
    double s=0;
    #pragma omp parallel for collapse(3) reduction(+:s) schedule(static)
    CELLS(a) s+=a(i,j,k);
    return d.sum(s)/d.count;
}
double remove_mean(Field& a,Domain& d) {
    double m=mean(a,d);
    #pragma omp parallel for collapse(3) schedule(static)
    CELLS(a) a(i,j,k)-=m;
    return m;
}
void apply_A(Field& x,Field& ax,Domain& d) {
    halo(x,d);
    double a=1/(d.h[0]*d.h[0]),b=1/(d.h[1]*d.h[1]),c=1/(d.h[2]*d.h[2]);
    #pragma omp parallel for collapse(3) schedule(static)
    CELLS(x) ax(i,j,k)=a*(2*x(i,j,k)-x(i-1,j,k)-x(i+1,j,k))+
                           b*(2*x(i,j,k)-x(i,j-1,k)-x(i,j+1,k))+
                           c*(2*x(i,j,k)-x(i,j,k-1)-x(i,j,k+1));
}
struct CGResult { int iterations=0; double rhs_rms=0,residual_rms=0,relative=0,removed_mean=0; };
struct Poisson {
    Field r,p,ap;
    Domain& d;
    const Options& o;
    int solves=0;
    Poisson(Domain& dom,const Options& opts):r(dom),p(dom),ap(dom),d(dom),o(opts) {}
    CGResult solve(Field& x,Field& b) {
        CGResult out;
        const auto solve_start=Clock::now();
        ++solves;
        out.removed_mean=remove_mean(b,d);
        remove_mean(x,d);
        const double b2=dot(b,b,d);
        out.rhs_rms=std::sqrt(b2/d.count);
        const double threshold=std::max(o.atol*std::sqrt(d.count),o.tol*std::sqrt(b2));
        auto true_residual=[&]() {
            apply_A(x,ap,d);
            #pragma omp parallel for collapse(3) schedule(static)
            CELLS(r) r(i,j,k)=b(i,j,k)-ap(i,j,k);
            return dot(r,r,d);
        };
        double rr=true_residual();
        remove_mean(r,d);
        #pragma omp parallel for collapse(3) schedule(static)
        CELLS(p) p(i,j,k)=r(i,j,k);
        rr=dot(r,r,d);
        bool converged=false;
        for(;;) {
            if(!std::isfinite(rr)) throw std::runtime_error("non-finite CG residual");
            if(std::sqrt(rr)<=threshold) {
                remove_mean(x,d);
                rr=true_residual();
                if(std::sqrt(rr)<=threshold) { converged=true; break; }
                // Restart from the actual residual if recursive residuals drift.
                remove_mean(r,d);
                #pragma omp parallel for collapse(3) schedule(static)
                CELLS(p) p(i,j,k)=r(i,j,k);
                rr=dot(r,r,d);
            }
            if(out.iterations>=o.max_iter) break;
            apply_A(p,ap,d);
            double pap=dot(p,ap,d);
            if(!std::isfinite(pap) || pap<=0) throw std::runtime_error("CG breakdown on zero-mean Neumann problem");
            double alpha=rr/pap;
            #pragma omp parallel for collapse(3) schedule(static)
            CELLS(x) { x(i,j,k)+=alpha*p(i,j,k); r(i,j,k)-=alpha*ap(i,j,k); }
            double next=dot(r,r,d), beta=next/rr;
            #pragma omp parallel for collapse(3) schedule(static)
            CELLS(p) p(i,j,k)=r(i,j,k)+beta*p(i,j,k);
            rr=next; ++out.iterations;
            // So that a run stopped by the wall time still gives a time per iteration.
            if(o.progress>0 && d.rank==0 && out.iterations%o.progress==0)
                std::cerr<<"PROGRESS solve="<<solves<<" iteration="<<out.iterations<<" solve_elapsed_s="<<elapsed(solve_start)
                         <<" residual_rms="<<std::sqrt(rr/d.count)<<std::endl;
        }
        out.residual_rms=std::sqrt(rr/d.count);
        out.relative=out.rhs_rms>0?out.residual_rms/out.rhs_rms:out.residual_rms;
        if(!converged) throw std::runtime_error("pressure CG did not converge within --max-iter");
        return out;
    }
};
using Velocity=std::array<Field,3>;
Velocity velocity(Domain& d) { return {Field(d),Field(d),Field(d)}; }
void velocity_halo(Velocity& u,Domain& d) { for(int a=0;a<3;++a) halo(u[a],d,a,true); }
double divergence(const Velocity& u,Field& div,Domain& d) {
    #pragma omp parallel for collapse(3) schedule(static)
    CELLS(div) div(i,j,k)=(u[0](i,j,k)-u[0](i-1,j,k))/d.h[0]+
                         (u[1](i,j,k)-u[1](i,j-1,k))/d.h[1]+
                         (u[2](i,j,k)-u[2](i,j,k-1))/d.h[2];
    return std::sqrt(dot(div,div,d)/d.count);
}
double maxabs(const Field& f,Domain& d) {
    double m=0; int invalid=0;
    #pragma omp parallel for collapse(3) reduction(max:m,invalid) schedule(static)
    CELLS(f) { double x=f(i,j,k); invalid=std::max(invalid,int(!std::isfinite(x))); m=std::max(m,std::abs(x)); }
    if(d.maximum(invalid)>0) throw std::runtime_error("non-finite field");
    return d.maximum(m);
}
void correct(Velocity& u,Field& pressure,Domain& d,double dt) {
    halo(pressure,d);
    for(int a=0;a<3;++a) {
        #pragma omp parallel for collapse(3) schedule(static)
        CELLS(pressure) {
            std::array<int,3> q{i,j,k},r=q; ++r[a];
            if(d.start[a]+q[a]==d.global[a]) u[a].at(q)=0;
            else u[a].at(q)-=dt*(pressure.at(r)-pressure.at(q))/d.h[a];
        }
    }
    velocity_halo(u,d);
}
void momentum(Velocity& u,Velocity& rhs,Domain& d,const Options& o,double time) {
    velocity_halo(u,d);
    for(int a=0;a<3;++a) {
        #pragma omp parallel for collapse(3) schedule(static)
        CELLS(u[a]) {
            std::array<int,3> q{i,j,k};
            if(d.start[a]+q[a]==d.global[a]) { rhs[a].at(q)=0; continue; }
            double adv=0,lap=0,qa=u[a].at(q);
            for(int b=0;b<3;++b) {
                auto plus=q,minus=q; ++plus[b]; --minus[b];
                double qp=u[a].at(plus),qm=u[a].at(minus),fp,fm;
                if(a==b) { fp=0.25*(qa+qp)*(qa+qp); fm=0.25*(qm+qa)*(qm+qa); }
                else {
                    auto cross=q,crossminus=minus; ++cross[a]; ++crossminus[a];
                    fp=0.25*(u[b].at(q)+u[b].at(cross))*(qa+qp);
                    fm=0.25*(u[b].at(minus)+u[b].at(crossminus))*(qm+qa);
                }
                adv+=(fp-fm)/d.h[b]; lap+=(qp-2*qa+qm)/(d.h[b]*d.h[b]);
            }
            double force=0;
            if(a==0) {
                double y=(d.start[1]+j-0.5)*d.h[1],z=(d.start[2]+k-0.5)*d.h[2];
                double kw=2*pi/o.length[1];
                double profile=o.profile=="uniform"?1:4*std::exp(-4*(1-z/o.length[2]))/(-std::expm1(-4.0));
                force=o.force/kw*(1+o.wind_alpha*std::sin(o.wind_omega*time))*std::cos(kw*y)*profile;
            }
            rhs[a].at(q)=-adv+o.nu*lap+force;
        }
    }
}
void initialize(Field& c,Domain& d,const Options& o) {
    #pragma omp parallel for collapse(3) schedule(static)
    CELLS(c) {
        double x=(d.start[0]+i-0.5)*d.h[0]/o.length[0];
        double y=(d.start[1]+j-0.5)*d.h[1]/o.length[1];
        double z=(d.start[2]+k-0.5)*d.h[2]/o.length[2];
        if(o.initial=="constant") c(i,j,k)=1;
        else if(o.initial=="cosine") c(i,j,k)=std::cos(pi*x)*std::cos(pi*y)*std::cos(pi*z);
        else c(i,j,k)=std::exp(-((x-.3)*(x-.3)+(y-.5)*(y-.5)+(z-.5)*(z-.5))/(2*.08*.08));
    }
}
void transport(Field& c,Field& next,const Velocity& u,Domain& d,const Options& o,double dt) {
    halo(c,d);
    #pragma omp parallel for collapse(3) schedule(static)
    CELLS(c) {
        std::array<int,3> q{i,j,k}; double value=c.at(q),change=0;
        for(int a=0;a<3;++a) {
            auto plus=q,minus=q; ++plus[a]; --minus[a];
            double up=u[a].at(q),um=u[a].at(minus),cp=c.at(plus),cm=c.at(minus);
            double fp=up*(up>=0?value:cp)-o.kappa*(cp-value)/d.h[a];
            double fm=um*(um>=0?cm:value)-o.kappa*(value-cm)/d.h[a];
            change+=(fp-fm)/d.h[a];
        }
        next.at(q)=value-dt*change;
    }
    c.swap(next);
}
double choose_dt(const Options& o,const Domain& d) {
    double a=0,s=0;
    for(double h:d.h) { a+=o.speed_bound/h; s+=1/(h*h); }
    return o.dt>0?o.dt:std::min({0.4/a,0.8/(a+2*o.kappa*s),o.nu>0?0.08/(o.nu*s):1e100});
}
void check_dt(const Velocity& u,Domain& d,const Options& o,double dt) {
    double a=0,s=0;
    for(int b=0;b<3;++b) { a+=maxabs(u[b],d)/d.h[b]; s+=1/(d.h[b]*d.h[b]); }
    if(dt*a>0.5 || dt*(a+2*o.kappa*s)>1 || o.nu*dt*s>0.1)
        throw std::runtime_error("fixed dt violates momentum/pigment stability bounds; reduce --dt / raise --speed-bound");
}
struct Stats { double mass=0,cmin=0,cmax=0,wmax=0,divmax=0,pmean=0; };
Stats stats(Field& c,Field& p,Velocity& u,Field& div,Domain& d) {
    Stats s; s.mass=mean(c,d)*d.count*d.h[0]*d.h[1]*d.h[2];
    double mn=std::numeric_limits<double>::infinity(),mx=-mn;
    #pragma omp parallel for collapse(3) reduction(min:mn) reduction(max:mx) schedule(static)
    CELLS(c) { mn=std::min(mn,c(i,j,k)); mx=std::max(mx,c(i,j,k)); }
    maxabs(c,d); s.cmin=-d.maximum(-mn); s.cmax=d.maximum(mx);
    s.wmax=maxabs(u[2],d); s.divmax=maxabs(div,d); s.pmean=mean(p,d);
    return s;
}
// Gather only for explicitly requested small validation/visualization outputs.
void write_field(const std::string& path,Field& c,Field& p,Velocity& u,Domain& d) {
    if(path.empty()) return;
    std::vector<double> local; local.reserve(static_cast<std::size_t>(c.nx)*c.ny*c.nz*8);
    CELLS(c) {
        for(int a=0;a<3;++a) local.push_back((d.start[a]+(a==0?i:(a==1?j:k))-.5)*d.h[a]);
        local.push_back(c(i,j,k)); local.push_back(p(i,j,k));
        for(int a=0;a<3;++a) local.push_back(u[a](i,j,k));
    }
    std::ofstream file;
    int bad=0;
    if(d.rank==0) { file.open(path); if(!file) bad=1; else file<<"x,y,z,c,phi,u_e,v_n,w_t\n"<<std::setprecision(17); }
    if(d.maximum(bad)>0) throw std::runtime_error("cannot open field output "+path);
    for(int rank=0;rank<d.ranks;++rank) {
        std::vector<double> buffer;
        if(d.rank==rank) buffer=local;
#ifdef USE_MPI
        if(local.size()>static_cast<std::size_t>(std::numeric_limits<int>::max())) throw std::runtime_error("field output exceeds MPI count");
        if(rank!=0) {
            int size=d.rank==rank?static_cast<int>(local.size()):0;
            if(d.rank==rank) { MPI_Send(&size,1,MPI_INT,0,90,d.cart); MPI_Send(local.data(),size,MPI_DOUBLE,0,91,d.cart); }
            if(d.rank==0) { MPI_Recv(&size,1,MPI_INT,rank,90,d.cart,MPI_STATUS_IGNORE); buffer.resize(size); MPI_Recv(buffer.data(),size,MPI_DOUBLE,rank,91,d.cart,MPI_STATUS_IGNORE); }
        }
#endif
        if(d.rank==0) for(std::size_t t=0;t<buffer.size();t+=8) {
            for(int k=0;k<8;++k) file<<(k?",":"")<<buffer[t+k];
            file<<'\n';
        }
    }
    if(d.rank==0) { file.flush(); bad=!file; }
    if(d.maximum(bad)>0) throw std::runtime_error("failed writing final field");
}

int execute(const Options& o) {
    Domain d(o);
    Field c(d),p(d),b(d),next(d),div(d);
    Velocity u=velocity(d),rhs=velocity(d),old=velocity(d),older=velocity(d);
    Poisson solver(d,o);
    double dt=choose_dt(o,d),mass0=0,loop_s=0,poisson_s=0,spinup_s=0,halo_s=0,reduction_s=0;
    double div_before=0,div_after=0,max_div=0,max_rel=0,error_rms=0,error_max=0;
    int total_iter=0,total_solves=0,max_iter=0;
    CGResult last;
    std::vector<std::string> traces;
    if(o.mode=="poisson-discrete" || o.mode=="poisson-continuous") {
        // Multiple Neumann eigenmodes avoid a trivial one-iteration-only check.
        #pragma omp parallel for collapse(3) schedule(static)
        CELLS(c) {
            double r[3]={(d.start[0]+i-.5)*d.h[0]/o.length[0],(d.start[1]+j-.5)*d.h[1]/o.length[1],(d.start[2]+k-.5)*d.h[2]/o.length[2]};
            double exact=0,forcing=0;
            const int modes[3][3]={{1,1,1},{2,1,0},{0,2,3}};
            const double weights[3]={1,.3,.2};
            for(int m=0;m<3;++m) {
                double v=weights[m],lambda=0;
                for(int a=0;a<3;++a) {
                    v*=std::cos(modes[m][a]*pi*r[a]);
                    double kh=modes[m][a]*pi/o.length[a];
                    lambda+=o.mode=="poisson-continuous"?kh*kh:4*std::pow(std::sin(kh*d.h[a]/2)/d.h[a],2);
                }
                exact+=v; forcing+=lambda*v;
            }
            c(i,j,k)=exact; b(i,j,k)=forcing;
        }
        d.barrier(); auto t=Clock::now(); last=solver.solve(p,b); poisson_s=loop_s=d.maximum(elapsed(t));
        double e2=0,em=0;
        #pragma omp parallel for collapse(3) reduction(+:e2) reduction(max:em) schedule(static)
        CELLS(c) { double e=p(i,j,k)-c(i,j,k); e2+=e*e; em=std::max(em,std::abs(e)); }
        error_rms=std::sqrt(d.sum(e2)/d.count); error_max=d.maximum(em);
        total_iter=last.iterations; total_solves=1; max_iter=total_iter; max_rel=last.relative;
    } else if(o.mode=="momentum-test") {
        // Analytic free-slip, divergence-free velocity with all three components.
        const double ka[3]={pi/o.length[0],pi/o.length[1],pi/o.length[2]};
        const double amplitude[3]={1,1,-(ka[0]+ka[1])/ka[2]};
        auto exact=[&](int component,const double* position) {
            double value=amplitude[component];
            for(int b=0;b<3;++b) value*=b==component?std::sin(ka[b]*position[b]):std::cos(ka[b]*position[b]);
            return value;
        };
        for(int a=0;a<3;++a) {
            #pragma omp parallel for collapse(3) schedule(static)
            CELLS(c) {
                int index[3]={i,j,k}; double position[3];
                for(int b=0;b<3;++b) position[b]=(d.start[b]+index[b]-(a==b?0:.5))*d.h[b];
                u[a](i,j,k)=d.start[a]+index[a]==d.global[a]?0:exact(a,position);
            }
        }
        Options test=o; test.force=0;
        momentum(u,rhs,d,test,0);
        double e2=0,em=0,count=0;
        for(int a=0;a<3;++a) {
            #pragma omp parallel for collapse(3) reduction(+:e2,count) reduction(max:em) schedule(static)
            CELLS(c) {
                int index[3]={i,j,k}; if(d.start[a]+index[a]==d.global[a]) continue;
                double position[3];
                for(int b=0;b<3;++b) position[b]=(d.start[b]+index[b]-(a==b?0:.5))*d.h[b];
                double adv=0,lap_factor=0;
                for(int b=0;b<3;++b) {
                    double derivative=amplitude[a]*ka[b];
                    for(int axis=0;axis<3;++axis) {
                        if(axis==b) derivative*=axis==a?std::cos(ka[axis]*position[axis]):-std::sin(ka[axis]*position[axis]);
                        else derivative*=axis==a?std::sin(ka[axis]*position[axis]):std::cos(ka[axis]*position[axis]);
                    }
                    adv+=exact(b,position)*derivative; lap_factor+=ka[b]*ka[b];
                }
                double expected=-adv-o.nu*lap_factor*exact(a,position);
                double e=rhs[a](i,j,k)-expected; e2+=e*e; em=std::max(em,std::abs(e)); ++count;
            }
        }
        error_rms=std::sqrt(d.sum(e2)/d.sum(count)); error_max=d.maximum(em);
        div_after=divergence(u,div,d); max_div=maxabs(div,d);
    } else if(o.mode=="projection") {
        // Deliberately divergent provisional velocity satisfying impermeable walls.
        for(int a=0;a<3;++a) {
            #pragma omp parallel for collapse(3) schedule(static)
            CELLS(c) {
                int q[3]={i,j,k}; double r[3];
                for(int t=0;t<3;++t) r[t]=(d.start[t]+q[t]-(t==a?0:.5))*d.h[t]/o.length[t];
                u[a](i,j,k)=d.start[a]+q[a]==d.global[a]?0:
                    (a+1)*.1*std::sin(pi*r[a])*(1+.3*std::cos(2*pi*r[(a+1)%3]))*(1+.2*std::cos(pi*r[(a+2)%3]));
            }
        }
        velocity_halo(u,d); div_before=divergence(u,div,d);
        #pragma omp parallel for collapse(3) schedule(static)
        CELLS(b) b(i,j,k)=-div(i,j,k)/dt;
        auto t=Clock::now(); last=solver.solve(p,b); poisson_s=d.maximum(elapsed(t));
        correct(u,p,d,dt); div_after=divergence(u,div,d); max_div=maxabs(div,d);
        // Check the exact discrete projection/residual identity, including removed RHS mean.
        apply_A(p,next,d);
        double err=0;
        #pragma omp parallel for collapse(3) reduction(max:err) schedule(static)
        CELLS(b) err=std::max(err,std::abs(div(i,j,k)+dt*(b(i,j,k)-next(i,j,k)+last.removed_mean)));
        error_max=d.maximum(err); total_iter=last.iterations; total_solves=1; max_rel=last.relative; max_iter=total_iter;
    } else {
        check_dt(u,d,o,dt);
        auto advance_flow=[&](int step) {
            momentum(u,rhs,d,o,step*dt);
            const double a0=step==0?1:(step==1?1.5:23.0/12);
            const double a1=step==0?0:(step==1?-.5:-16.0/12),a2=step<2?0:5.0/12;
            for(int a=0;a<3;++a) {
                #pragma omp parallel for collapse(3) schedule(static)
                CELLS(c) u[a](i,j,k)+=dt*(a0*rhs[a](i,j,k)+a1*old[a](i,j,k)+a2*older[a](i,j,k));
                older[a].swap(old[a]); old[a].swap(rhs[a]);
            }
            velocity_halo(u,d); div_before=divergence(u,div,d);
            #pragma omp parallel for collapse(3) schedule(static)
            CELLS(b) b(i,j,k)=-div(i,j,k)/dt;
            auto t=Clock::now(); last=solver.solve(p,b); poisson_s+=elapsed(t);
            correct(u,p,d,dt); div_after=divergence(u,div,d);
            max_div=std::max(max_div,maxabs(div,d)); max_rel=std::max(max_rel,last.relative);
            total_iter+=last.iterations; ++total_solves; max_iter=std::max(max_iter,last.iterations);
            if(div_after>10*(dt*(std::max(o.atol,o.tol*last.rhs_rms)+std::abs(last.removed_mean))+1e-13))
                throw std::runtime_error("projection failed divergence bound");
            check_dt(u,d,o,dt);
        };
        d.barrier(); auto sp=Clock::now();
        for(int n=0;n<o.spinup;++n) advance_flow(n);
        spinup_s=d.maximum(elapsed(sp));
        initialize(c,d,o); mass0=mean(c,d)*d.count*d.h[0]*d.h[1]*d.h[2];
        total_iter=total_solves=max_iter=0; max_div=max_rel=poisson_s=0;
        d.comm_s=d.reduce_s=0; d.barrier(); auto t=Clock::now();
        for(int n=0;n<o.steps;++n) {
            transport(c,next,u,d,o,dt); advance_flow(o.spinup+n);
            if(!o.trace.empty()) {
                auto s=stats(c,p,u,div,d);
                if(d.rank==0) {
                    std::ostringstream line; line<<std::setprecision(17)<<n+1<<','<<(n+1)*dt<<','<<(o.spinup+n+1)*dt<<','
                        <<last.iterations<<','<<last.relative<<','<<last.residual_rms<<','<<last.removed_mean<<','
                        <<div_before<<','<<div_after<<','<<s.divmax<<','<<s.mass<<','<<s.cmin<<','<<s.cmax<<','<<s.wmax;
                    traces.push_back(line.str());
                }
            }
        }
        const double local_loop_s=elapsed(t),local_halo_s=d.comm_s,local_reduce_s=d.reduce_s;
        loop_s=d.maximum(local_loop_s); poisson_s=d.maximum(poisson_s);
        halo_s=d.maximum(local_halo_s); reduction_s=d.maximum(local_reduce_s);
    }
    if(o.mode!="simulate") {
        const double local_halo_s=d.comm_s,local_reduce_s=d.reduce_s;
        halo_s=d.maximum(local_halo_s); reduction_s=d.maximum(local_reduce_s);
    }
    auto s=stats(c,p,u,div,d);
    double pigment_error=0;
    if(o.mode=="simulate" && o.force==0 && (o.initial=="cosine" || o.initial=="constant")) {
        initialize(next,d,o); double lambda=0;
        for(int a=0;a<3;++a) lambda+=4*o.kappa*std::pow(std::sin(pi/(2*o.n[a]))/d.h[a],2);
        const double factor=o.initial=="cosine"?std::pow(1-dt*lambda,o.steps):1;
        double e=0;
        #pragma omp parallel for collapse(3) reduction(max:e) schedule(static)
        CELLS(c) e=std::max(e,std::abs(c(i,j,k)-factor*next(i,j,k)));
        pigment_error=d.maximum(e);
    }
    auto io=Clock::now(); write_field(o.output,c,p,u,d);
    if(d.rank==0 && !o.trace.empty()) {
        std::ofstream f(o.trace);
        f<<"step,pigment_time,flow_time,cg_iterations,relative_residual,residual_rms,rhs_mean_removed,div_before_rms,div_after_rms,div_after_max,mass,c_min,c_max,w_max\n";
        for(const auto& row:traces) f<<row<<'\n';
        f.flush(); if(!f) throw std::runtime_error("cannot write trace "+o.trace);
    }
    double output_s=d.maximum(elapsed(io));
    int threads=1;
#ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp single
        threads=omp_get_num_threads();
    }
#endif
    if(d.rank==0) std::cout<<std::setprecision(17)<<"RESULT {\"model\":\"B-pressure-projection\",\"mode\":\""<<o.mode
        <<"\",\"ranks\":"<<d.ranks<<",\"threads\":"<<threads<<",\"nodes\":"<<d.nodes
        <<",\"nx\":"<<o.n[0]<<",\"ny\":"<<o.n[1]<<",\"nz\":"<<o.n[2]
        <<",\"px\":"<<d.dims[0]<<",\"py\":"<<d.dims[1]<<",\"pz\":"<<d.dims[2]
        <<",\"pressure_unknowns\":"<<std::fixed<<std::setprecision(0)<<d.count<<std::defaultfloat<<std::setprecision(17)
        <<",\"dt\":"<<dt<<",\"steps\":"<<o.steps<<",\"spinup\":"<<o.spinup
        <<",\"loop_s\":"<<loop_s<<",\"poisson_s\":"<<poisson_s<<",\"spinup_s\":"<<spinup_s<<",\"output_s\":"<<output_s
        <<",\"halo_s\":"<<halo_s<<",\"reduction_s\":"<<reduction_s
        <<",\"cg_solves\":"<<total_solves<<",\"cg_iterations\":"<<total_iter<<",\"cg_max_iterations\":"<<max_iter
        <<",\"max_relative_residual\":"<<max_rel<<",\"last_residual_rms\":"<<last.residual_rms
        <<",\"rhs_mean_removed\":"<<last.removed_mean<<",\"div_before_rms\":"<<div_before
        <<",\"div_after_rms\":"<<div_after<<",\"max_divergence\":"<<max_div
        <<",\"pressure_mean\":"<<s.pmean<<",\"mass_initial\":"<<mass0<<",\"mass_final\":"<<s.mass
        <<",\"mass_abs_error\":"<<std::abs(s.mass-mass0)<<",\"c_min\":"<<s.cmin<<",\"c_max\":"<<s.cmax
        <<",\"w_max\":"<<s.wmax<<",\"error_rms\":"<<error_rms<<",\"error_max\":"<<error_max
        <<",\"pigment_discrete_error\":"<<pigment_error<<"}\n";
    return 0;
}
int main(int argc,char** argv) {
#ifdef USE_MPI
    int provided=0; MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
    if(provided<MPI_THREAD_FUNNELED) MPI_Abort(MPI_COMM_WORLD,2);
#endif
    int status=0;
    try { auto o=parse(argc,argv); if(o.steps!=-1) status=execute(o); }
    catch(const std::exception& e) {
        std::cerr<<"ERROR: "<<e.what()<<'\n'; status=1;
#ifdef USE_MPI
        MPI_Abort(MPI_COMM_WORLD,1);
#endif
    }
#ifdef USE_MPI
    MPI_Finalize();
#endif
    return status;
}
