#!/bin/bash
module purge
module load pgi
module load cuda/7.5

export MPI_HOME=$HOME/APPS/OPENMPI/1.10.2/pgi-16.5
export PATH=$MPI_HOME/bin:$PATH
export LD_LIBRARY_PATH=$MPI_HOME/lib:$LD_LIBRARY_PATH

export HWLOC_HOME=$HOME/APPS/HWLOC/1.11.4
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HWLOC_HOME/lib
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$HWLOC_HOME/lib/pkgconfig
export PATH=$PATH:$HWLOC_HOME/bin
export MANPATH=$MANPATH:$HWLOC_HOME/share/man
