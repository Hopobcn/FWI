#!/bin/bash

###############################################################
########### DO NOT MODIFY WALL-CLOCK-TIME LIMIT! ##############
# @ job_name         = fwi.GNU.seq
# @ initialdir       = .
# @ output           = fwi.GNU.omp.%j.out
# @ error            = fwi.GNU.omp.%j.err
# @ total_tasks      = 1
# @ tasks_per_node   = 1
# @ cpus_per_task    = 16
# @ gpus_per_node    = 4
# @ wall_clock_limit = 00:05:00
# @ features         = k80
###############################################################
ulimit -c unlimited


source ../environment_gcc.sh

# OpenMP num threads and thread affinity:
export OMP_NUM_THREADS=16

# For gcc that support OPenMP 3.1 or newer use this ones:
export OMP_PROC_BIND=true
# otherwise use this one:
#export GOMP_CPU_AFFINITY="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15"

# Estimated elapsed time: 3:00 minutes
../../build/bin/fwi ../../data/fwi_params.txt ../../data/fwi_frequencies.medium.txt

mv fwi.00.log fwi.GNU.omp.${SLURM_JOB_ID}.log
