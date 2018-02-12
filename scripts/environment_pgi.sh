#!/bin/bash
module purge
module load pgi/17.10
module load gcc/5.1.0    # CUDA 8.0 requires gcc 5 or previous
module load cuda/8.0     # CUDA 9.0 requires gcc 6 or previous
module load cmake/3.9.6

# Make CTEST VERBOSE when some test fails:
export CTEST_OUTPUT_ON_FAILURE=1
