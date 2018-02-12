#!/bin/bash
module purge
module load intel
module load cmake/3.9.6

# Make CTEST VERBOSE when some test fails:
export CTEST_OUTPUT_ON_FAILURE=1
