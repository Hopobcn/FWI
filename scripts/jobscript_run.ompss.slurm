#!/bin/bash

##
## Run this script using 'make run-openmp' from your build directory
## (if you want to use this script manually paths should be modified)
##

###############################################################
########### DO NOT MODIFY WALL-CLOCK-TIME LIMIT! ##############
#SBATCH -J fwi-ompss
#SBATCH -D .
#SBATCH --output=fwi.ompss.%j.out
#SBATCH --error=fwi.ompss.%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --constraint=k80
#SBATCH --time=01:20:00
###############################################################

ulimit -c unlimited

echo "PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}" # FWI project directory
echo "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}" # directory where 'fwi' binary is
echo "COMPILER_ID:        ${COMPILER_ID}"        # compiler used

if [ "$COMPILER_ID" == "GNU" ]
then
    source ${PROJECT_SOURCE_DIR}/scripts/environment_mcc.sh

    # OmpSs num threads and thread affinity:
    export NX_ARGS='--smp-workers=16 --cpus-per-socket=8 --num-sockets=2 --summary'
    #export OMP_NUM_THREADS=16
    # For gcc that support OPenMP 3.1 or newer use this ones:
    #export OMP_PROC_BIND=true
    # otherwise use this one:
    #export GOMP_CPU_AFFINITY="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15"

    export NX_INSTRUMENTATION=extrae
    export EXTRAE_CONFIG_FILE=${PROJECT_SOURCE_DIR}/scripts/extrae.xml

elif [ "$COMPILER_ID" == "Intel" ]
then
    source ${PROJECT_SOURCE_DIR}/scripts/environment_imcc.sh

    # OpenMP num threads and thread affinity:
    #export OMP_NUM_THREADS=16
    #export KMP_AFFINITY=verbose,compact
    export NX_ARGS='--smp-workers=16 --cpus-per-socket=8 --num-sockets=2 --summary'

    export NX_INSTRUMENTATION=extrae
    export EXTRAE_CONFIG_FILE=${PROJECT_SOURCE_DIR}/scripts/extrae.xml

else
    echo "ERROR: COMPILER_ID ${COMPILER_ID} not supported"
fi

${PROJECT_BINARY_DIR}/fwi ${PROJECT_SOURCE_DIR}/data/fwi_params.txt ${PROJECT_SOURCE_DIR}/data/fwi_frequencies.txt

mv fwi.00.log fwi.ompss.${SLURM_JOB_ID}.log
