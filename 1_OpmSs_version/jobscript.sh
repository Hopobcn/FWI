#!/bin/bash
##16+1
#BSUB -n 1
#BSUB -oo output.out
#BSUB -eo output.err
#BSUB -R "span[ptile=16]"
#BSUB -x
#BSUB -W 06:00
#BSUB -q bsc_cs
#BSUB -J fwi.ompss
ulimit -c unlimited

date
source environment.sh

# Intentar imitar /gpfs/scratch/bsc15/bsc15685/offload/apps/offload.test/mn3/mn3_x2m.lsf
date

export OMP_NUM_THREADS=16
export NX_ARGS="--instrument-default=omp --summary --enable-yield"

date

./trace.sh ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

date
