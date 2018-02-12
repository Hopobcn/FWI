#!/bin/bash

##
## Run this script using 'make run-mpi+openacc' from your build directory
## (if you want to use this script manually paths should be modified)
##

###############################################################
########### DO NOT MODIFY WALL-CLOCK-TIME LIMIT! ##############
#SBATCH -J fwi-mpi+acc
#SBATCH -D .
#SBATCH --output=fwi.mpi+acc.%j.out
#SBATCH --error=fwi.mpi+acc.%j.err
#SBATCH --ntasks=2
#SBATCH --gres gpu:2
#SBATCH --cpus-per-task=4
#SBATCH --constraint=k80
#SBATCH --time=00:05:00
###############################################################

ulimit -c unlimited

echo "PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}" # FWI project directory
echo "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}" # directory where 'fwi' binary is
echo "COMPILER_ID:        ${COMPILER_ID}"        # compiler used

if [ "$COMPILER_ID" == "PGI" ]
then
    source ${PROJECT_SOURCE_DIR}/scripts/environment_pgi.sh

elif [ "$COMPILER_ID" == "GNU" ]
then
    echo "ERROR: GCC does not have OpenACC 2.5 support"

elif [ "$COMPILER_ID" == "Intel" ]
then
    echo "ERROR: Intel does not have OpenACC support"
else
    echo "ERROR: COMPILER_ID not recognized as PGI/GNU/Intel"
fi

# In case you compiled with openmp, force only one thread
# If your MULTI-GPU implementation uses OpenMP, modify this!
export OMP_NUM_THREADS=1

export FWIDIR=${PROJECT_SOURCE_DIR}

mpirun -np 2 ${PROJECT_BINARY_DIR}/fwi fwi_schedule.txt

mv fwi.00.log fwi.PGI.mpi.${SLURM_JOB_ID}.00.log
mv fwi.01.log fwi.PGI.mpi.${SLURM_JOB_ID}.01.log
