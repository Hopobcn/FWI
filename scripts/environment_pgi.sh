#!/bin/bash
module purge
module load pgi/16.9
module load cuda/7.5
#module load openmpi/1.10.2_cuda_pgi
module load gcc/4.9.3     # NVCC requires a GCC version that understands -std=c++11
module load cmake/3.6.2

# Make CTEST VERBOSE when some test fails:
export CTEST_OUTPUT_ON_FAILURE=1
