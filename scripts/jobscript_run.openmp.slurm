#!/bin/bash

##
## Run this script using 'make run-openmp' from your build directory
## (if you want to use this script manually paths should be modified)
##

###############################################################
########### DO NOT MODIFY WALL-CLOCK-TIME LIMIT! ##############
#SBATCH -J fwi-omp
#SBATCH -D .
#SBATCH --output=fwi.omp.%j.out
#SBATCH --error=fwi.omp.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --constraint=k80
#SBATCH --time=00:20:00
###############################################################

ulimit -c unlimited

echo "PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}" # FWI project directory
echo "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}" # directory where 'fwi' binary is
echo "COMPILER_ID:        ${COMPILER_ID}"        # compiler used

if [ "$COMPILER_ID" == "PGI" ]
then
    source ${PROJECT_SOURCE_DIR}/scripts/environment_pgi.sh

    # OpenMP num threads and thread affinity:
    export OMP_NUM_THREADS=16
    export MP_BIND=yes
    export MP_BLIST=0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15

elif [ "$COMPILER_ID" == "GNU" ]
then
    source ${PROJECT_SOURCE_DIR}/scripts/environment_gcc.sh

    # OpenMP num threads and thread affinity:
    export OMP_NUM_THREADS=16
    # For gcc that support OPenMP 3.1 or newer use this ones:
    export OMP_PROC_BIND=true
    # otherwise use this one:
    #export GOMP_CPU_AFFINITY="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15"

elif [ "$COMPILER_ID" == "Intel" ]
then
    source ${PROJECT_SOURCE_DIR}/scripts/environment_icc.sh

    # OpenMP num threads and thread affinity:
    export OMP_NUM_THREADS=16
    export KMP_AFFINITY=verbose,compact

else
    echo "ERROR: COMPILER_ID not recognized as PGI/GNU/Intel"
fi

export FWIDIR=${PROJECT_SOURCE_DIR}

${PROJECT_BINARY_DIR}/fwi fwi_schedule.txt

mv fwi.00.log fwi.omp.${SLURM_JOB_ID}.log
