#!/bin/bash

###############################################################
########### DO NOT MODIFY WALL-CLOCK-TIME LIMIT! ##############
# @ job_name         = fwi.PGI.seq
# @ initialdir       = .
# @ output           = fwi.PGI.omp.%j.out
# @ error            = fwi.PGI.omp.%j.err
# @ total_tasks      = 1
# @ tasks_per_node   = 1
# @ cpus_per_task    = 16
# @ gpus_per_node    = 4
# @ wall_clock_limit = 00:20:00
# @ features         = k80
###############################################################

ulimit -c unlimited

source ../environment_pgi.sh

# OpenMP num threads and thread affinity:
export OMP_NUM_THREADS=16
export MP_BIND=yes
export MP_BLIST=0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15

# estimated elapsed time 17:00 minutes
../../build/bin/fwi ../../data/fwi_params.txt ../../data/fwi_frequencies.txt

mv fwi.00.log fwi.PGI.omp.${SLURM_JOB_ID}.log
