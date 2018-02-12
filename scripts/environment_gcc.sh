#!/bin/bash
module purge
module load gcc/7.1.0
module load cmake/3.9.6

# Make CTEST VERBOSE when some test fails:
export CTEST_OUTPUT_ON_FAILURE=1
