#!/bin/bash
module purge
module load intel
module load impi
module load cmake/3.6.2

# Make CTEST VERBOSE when some test fails:
export CTEST_OUTPUT_ON_FAILURE=1
