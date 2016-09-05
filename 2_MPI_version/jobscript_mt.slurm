#!/bin/bash

# @ job_name      = fwi.mpi
# @ initialdir    = .
# @ output        = Results/%j.out
# @ error         = Results/%j.err
# @ total_tasks   = 8
# @ cpus_per_task = 2
# @ gpus_per_node = 4
# @ wall_clock_limit = 10:00:00
# @ features = k80
ulimit -c unlimited

NPROCS=16
source environment_mt.sh
export OMP_NUM_THREADS=2

mpirun -np 8 ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt
