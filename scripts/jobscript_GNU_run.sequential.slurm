#!/bin/bash

###############################################################
########### DO NOT MODIFY WALL-CLOCK-TIME LIMIT! ##############
# @ job_name         = fwi.GNU.seq
# @ initialdir       = .
# @ output           = fwi.GNU.seq.%j.out
# @ error            = fwi.GNU.seq.%j.err
# @ total_tasks      = 1
# @ tasks_per_node   = 1
# @ cpus_per_task    = 1
# @ gpus_per_node    = 1
# @ wall_clock_limit = 00:08:00
# @ features         = k80
###############################################################

ulimit -c unlimited

source ../environment_gcc.sh

# In case you compiled with openmp, force only one thread
export OMP_NUM_THREADS=1

# Estimated elapsed time: 6:00 minutes
../../build/bin/fwi ../../data/fwi_params.txt ../../data/fwi_frequencies.profile.txt

mv fwi.00.log fwi.GNU.seq.${SLURM_JOB_ID}.log
