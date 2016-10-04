#!/bin/bash

###############################################################
########### DO NOT MODIFY WALL-CLOCK-TIME LIMIT! ##############
# @ job_name         = fwi.Intel.omp
# @ initialdir       = .
# @ output           = fwi.Intel.omp.%j.out
# @ error            = fwi.Intel.omp.%j.err
# @ total_tasks      = 1
# @ tasks_per_node   = 1
# @ cpus_per_task    = 16
# @ gpus_per_node    = 4
# @ wall_clock_limit = 00:05:00
# @ features         = k80
###############################################################

ulimit -c unlimited

# This script should be called from the FWI project directory

source ../environment_icc.sh

# OpenMP num threads and thread affinity:
export OMP_NUM_THREADS=16
export KMP_AFFINITY=verbose,compact

# expected elapsed time: 2:30 minutes
../../build/bin/fwi ../../data/fwi_params.txt ../../data/fwi_frequencies.medium.txt

mv fwi.00.log fwi.Intel.omp.${SLURM_JOB_ID}.log
