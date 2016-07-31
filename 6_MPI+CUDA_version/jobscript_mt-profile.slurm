#!/bin/bash

# @ job_name      = fwi.mpi-openacc-cuda
# @ initialdir    = .
# @ output        = Results/%j.out
# @ error         = Results/%j.err
# @ total_tasks   = 4
# @ cpus_per_task = 4
# @ gpus_per_node = 4
# @ wall_clock_limit = 00:10:00
# @ features = k80
ulimit -c unlimited

NPROCS=4
source environment_mt.sh
export OMP_NUM_THREADS=4

mpirun -np ${NPROCS} nvprof -o fwi.%q{OMPI_COMM_WORLD_RANK}.nvprof ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

