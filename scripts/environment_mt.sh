#!/bin/bash
module purge
module load pgi
module load cuda/7.5
module load openmpi/1.10.2_cuda_pgi
module load cmake/3.6.2

export HWLOC_HOME=$HOME/APPS/HWLOC/1.11.4
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HWLOC_HOME/lib
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$HWLOC_HOME/lib/pkgconfig
export PATH=$PATH:$HWLOC_HOME/bin
export MANPATH=$MANPATH:$HWLOC_HOME/share/man

# Make CTEST VERBOSE when some test fails:
export CTEST_OUTPUT_ON_FAILURE=1
