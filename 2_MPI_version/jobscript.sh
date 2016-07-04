#!/bin/bash
#BSUB -n 2
#BSUB -oo fwi_mpi.out
#BSUB -eo fwi_mpi.err
#BSUB -R "span[ptile=16]"
#BSUB -x
#BSUB -W 00:10
#BSUB -q debug
#BSUB -J fwi.mpi
ulimit -c unlimited

date
source environment.sh

# Intentar imitar /gpfs/scratch/bsc15/bsc15685/offload/apps/offload.test/mn3/mn3_x2m.lsf
date

export OMP_NUM_THREADS=16

date

mpirun -n 2 ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

date
