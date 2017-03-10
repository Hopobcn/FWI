#!/bin/bash
module purge
module load gcc/6.1.0
module load openmpi/1.10.2
module load cmake/3.6.2
module load extrae
module load ompss

# extrae requires papi
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/bullxde/perftools/papi-bull/5.4.0.0/lib64

# Make CTEST VERBOSE when some test fails:
export CTEST_OUTPUT_ON_FAILURE=1
