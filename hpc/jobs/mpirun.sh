#!/bin/bash
# mpirun for Vanda's Open MPI 5: it starts its daemons over ssh, which needs
# the PRRTE install's path and the modules' libraries (hpc/README.md, "Experiment jobs").
exec mpirun --prefix "$EBROOTPRRTE" -x LD_LIBRARY_PATH "$@"
