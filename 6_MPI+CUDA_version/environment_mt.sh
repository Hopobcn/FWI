#!/bin/bash
module purge
module load pgi
module load cuda/7.5

export MPI_HOME=$HOME/APPS/OPENMPI/1.10.2/pgi-16.5
export PATH=$MPI_HOME/bin:$PATH
export LD_LIBRARY_PATH=$MPI_HOME/lib:$LD_LIBRARY_PATH
#
