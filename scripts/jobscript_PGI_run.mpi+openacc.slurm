#!/bin/bash

###############################################################
########### DO NOT MODIFY WALL-CLOCK-TIME LIMIT! ##############
# @ job_name         = fwi.PGI.mpi+openacc
# @ initialdir       = .
# @ output           = fwi.PGI.mpi+acc.%j.out
# @ error            = fwi.PGI.mpi+acc.%j.err
# @ total_tasks      = 2
# @ tasks_per_node   = 2
# @ cpus_per_task    = 4
# @ gpus_per_node    = 2
# @ wall_clock_limit = 00:05:00
# @ features         = k80
###############################################################

ulimit -c unlimited

source ../environment_pgi.sh

mpirun -np 2 ../../build/bin/fwi ../../data/fwi_params.txt ../../data/fwi_frequencies.medium.txt

mv fwi.00.log fwi.PGI.mpi.${SLURM_JOB_ID}.00.log
mv fwi.01.log fwi.PGI.mpi.${SLURM_JOB_ID}.01.log
