#!/bin/bash

# @ job_name      = fwi.openacc
# @ initialdir    = .
# @ output        = output.out
# @ error         = output.err
# @ total_tasks   = 1
# @ cpus_per_task = 6
# @ gpus_per_node = 1
# @ wall_clock_limit = 00:00:10
# @ features = k80
ulimit -c unlimited

source environment_mt.sh
export OMP_NUM_THREADS=16

./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

