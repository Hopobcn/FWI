#!/bin/bash
module purge
module load gcc/6.1.0
module load openmpi/1.10.2
module load cmake/3.6.2

# Make CTEST VERBOSE when some test fails:
export CTEST_OUTPUT_ON_FAILURE=1
