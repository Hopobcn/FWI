#!/bin/bash

###############################################################
########### DO NOT MODIFY WALL-CLOCK-TIME LIMIT! ##############
# @ job_name         = fwi.PGI.acc
# @ initialdir       = .
# @ output           = fwi.PGI.acc.%j.out
# @ error            = fwi.PGI.acc.%j.err
# @ total_tasks      = 1
# @ tasks_per_node   = 1
# @ cpus_per_task    = 4
# @ gpus_per_node    = 1
# @ wall_clock_limit = 00:05:00
# @ features         = k80
###############################################################

ulimit -c unlimited

source ../environment_pgi.sh

# estimated elapsed time 85 seconds
../../build/bin/fwi ../../data/fwi_params.txt ../../data/fwi_frequencies.txt

mv fwi.00.log fwi.PGI.acc.${SLURM_JOB_ID}.log
