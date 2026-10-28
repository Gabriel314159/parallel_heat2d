# Technical references

These sources informed the planned design. Record any additional algorithmic or implementation references used during development.

- ENCCS, [Heat diffusion mini-app](https://enccs.github.io/openmp-gpu/miniapp/): five-point finite differences, explicit time stepping and field storage.
- ENCCS, [Non-blocking point-to-point](https://enccs.github.io/intermediate-mpi/non-blocking-communication-pt1/): halo exchange, nonblocking messages and communication/computation overlap.
- ENCCS, [MPI and threads in practice](https://enccs.github.io/intermediate-mpi/mpi-and-threads-pt2/): hybrid MPI/OpenMP organization and thread placement.
- MPI Forum, [MPI 4.1 standard](https://www.mpi-forum.org/docs/mpi-4.1/mpi41-report/mpi41-report.htm): message passing, thread support, topology and timing semantics.
- CMake, [FindMPI](https://cmake.org/cmake/help/latest/module/FindMPI.html) and [FindOpenMP](https://cmake.org/cmake/help/latest/module/FindOpenMP.html): imported build targets and dependency discovery.

Numerical method and physics:

- R. J. LeVeque, *Finite Volume Methods for Hyperbolic Problems*, Cambridge University Press, 2002: conservative flux form and upwind advection.
- G. K. Vallis, *Atmospheric and Oceanic Fluid Dynamics*, 2nd ed., Cambridge University Press, 2017: the vorticity–streamfunction form of two-dimensional flow and wind-driven gyres.
- A. Arakawa, "Computational design for long-term numerical integration of the equations of fluid motion: two-dimensional incompressible flow. Part I", *Journal of Computational Physics* 1 (1966) 119–143: the energy- and enstrophy-conserving Jacobian.
- J. R. Shewchuk, *An Introduction to the Conjugate Gradient Method Without the Agonizing Pain*, Carnegie Mellon University, 1994: the conjugate gradient method used for the Poisson equation.
- Y. Saad, *Iterative Methods for Sparse Linear Systems*, 2nd ed., SIAM, 2003: convergence of CG and preconditioning.

The source code will use standard MPI C calls from C++, not the removed legacy MPI C++ bindings. Platform module names, launch commands and queue settings must come from current Atlas/Vanda documentation.
