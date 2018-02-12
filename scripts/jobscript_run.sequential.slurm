#!/bin/bash

##
## Run this script using 'make run-seq' from your build directory
## (if you want to use this script manually paths should be modified)
##

###############################################################
########### DO NOT MODIFY WALL-CLOCK-TIME LIMIT! ##############
#SBATCH -J fwi-seq
#SBATCH -D .
#SBATCH --output=fwi.seq.%j.out
#SBATCH --error=fwi.seq.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --constraint=k80
#SBATCH --time=00:10:00
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
    source ${PROJECT_SOURCE_DIR}/scripts/environment_gcc.sh
elif [ "$COMPILER_ID" == "Intel" ]
then
    source ${PROJECT_SOURCE_DIR}/scripts/environment_icc.sh
else
    echo "ERROR: COMPILER_ID not recognized as PGI/GNU/Intel"
fi

# In case you compiled with openmp, force only one thread
export OMP_NUM_THREADS=1

export FWIDIR=${PROJECT_SOURCE_DIR}


${PROJECT_BINARY_DIR}/fwi fwi_schedule.txt

mv fwi.00.log fwi.seq.${SLURM_JOB_ID}.log
