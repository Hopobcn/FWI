#!/bin/bash

##
## Run this script using 'make irun' from your build directory
##

ulimit -c unlimited

echo "PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}" # FWI project directory
echo "PROJECT_BINARY_DIR: ${PROJECT_BINARY_DIR}" # directory where 'fwi' binary is
echo "COMPILER_ID:        ${COMPILER_ID}"        # compiler used
echo "---"
echo "${PROJECT_BINARY_DIR}/fwi fwi_schedule.txt"
echo "---"

export FWIDIR=${PROJECT_SOURCE_DIR}
${PROJECT_BINARY_DIR}/fwi fwi_schedule.txt
