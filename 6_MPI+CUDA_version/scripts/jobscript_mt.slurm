#!/bin/bash

# @ job_name      = fwi.mpi-openacc-cuda
# @ initialdir    = .
# @ output        = Results/%j.out
# @ error         = Results/%j.err
# @ total_tasks   = 4
# @ tasks_per_node= 4
# @ cpus_per_task = 1
# @ gpus_per_node = 4
# @ wall_clock_limit = 01:00:00
# @ features = k80
ulimit -c unlimited

NTASKS=4
source environment_mt.sh
export OMP_NUM_THREADS=1

mpirun -np $NTASKS ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

