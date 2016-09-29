#!/bin/bash
module purge
module load pgi/16.5
module load cuda/7.5
module load openmpi/1.10.2_cuda_pgi
module load cmake/3.6.2

# Make CTEST VERBOSE when some test fails:
export CTEST_OUTPUT_ON_FAILURE=1
