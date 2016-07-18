#!/bin/bash

# @ job_name      = fwi.openacc-cuda
# @ initialdir    = .
# @ output        = 1/%j.out
# @ error         = 1/%j.err
# @ total_tasks   = 1
# @ cpus_per_task = 16
# @ gpus_per_node = 4
# @ wall_clock_limit = 00:01:00
# @ features = k80
ulimit -c unlimited

source environment_mt.sh
export OMP_NUM_THREADS=4

./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

