#!/bin/bash
# mpirun for Vanda's Open MPI 5: it starts its daemons over ssh, which needs
# the PRRTE install's path (hpc/README.md, "Experiment jobs").
exec mpirun --prefix "$EBROOTPRRTE" "$@"
