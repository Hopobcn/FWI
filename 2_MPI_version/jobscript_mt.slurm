#!/bin/bash

# @ job_name      = fwi.mpi
# @ initialdir    = .
# @ output        = Results/%j.out
# @ error         = Results/%j.err
# @ total_tasks   = 16
# @ cpus_per_task = 1
# @ gpus_per_node = 4
# @ wall_clock_limit = 17:00:00
# @ features = k80
ulimit -c unlimited

NPROCS=16
source environment_mt.sh
export OMP_NUM_THREADS=1

#mpirun -np 16 ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt
mpirun -np 16 ./fwi.intel64 ../SetupParams/fwi_params.txt ./fwi_frequencies.4.txt
