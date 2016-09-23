#!/bin/bash

# @ job_name      = fwi.openacc-cuda
# @ initialdir    = .
# @ output        = Results/%j.out
# @ error         = Results/%j.err
# @ total_tasks   = 1
# @ cpus_per_task = 4
# @ gpus_per_node = 4
# @ wall_clock_limit = 05:00:00
# @ features = k80
ulimit -c unlimited

source environment_mt.sh
export OMP_NUM_THREADS=4

./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

