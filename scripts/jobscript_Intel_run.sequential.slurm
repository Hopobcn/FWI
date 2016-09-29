#!/bin/bash

###############################################################
########### DO NOT MODIFY WALL-CLOCK-TIME LIMIT! ##############
# @ job_name         = fwi.Intel.seq
# @ initialdir       = .
# @ output           = fwi.Intel.seq.%j.out
# @ error            = fwi.Intel.seq.%j.err
# @ total_tasks      = 1
# @ tasks_per_node   = 1
# @ cpus_per_task    = 1
# @ gpus_per_node    = 1
# @ wall_clock_limit = 00:05:00
# @ features         = k80
###############################################################

ulimit -c unlimited

source ../environment_icc.sh

# In case you compiled with openmp, force only one thread
export OMP_NUM_THREADS=1

../../build/bin/fwi ../../data/fwi_params.txt ../../data/fwi_frequencies.small.txt
