add_executable(fwi
    fwi_main.c
)

target_link_libraries(fwi
    fwi-core
    m
)

add_custom_target(irun
    COMMAND PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR} PROJECT_BINARY_DIR=${CMAKE_RUNTIME_OUTPUT_DIRECTORY} COMPILER_ID=${CMAKE_C_COMPILER_ID} ${PROJECT_SOURCE_DIR}/scripts/irun.sh
    DEPENDS fwi
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/scripts/output
    COMMENT "outputs will be in ${PROJECT_SOURCE_DIR}/scripts/output/"
    VERBATIM
    USES_TERMINAL
)

add_custom_target(run-seq
    COMMAND sbatch --export=PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR},PROJECT_BINARY_DIR=${CMAKE_RUNTIME_OUTPUT_DIRECTORY},COMPILER_ID=${CMAKE_C_COMPILER_ID} ${PROJECT_SOURCE_DIR}/scripts/jobscript_run.sequential.slurm
    DEPENDS fwi
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/scripts/output
    COMMENT "outputs will be in ${PROJECT_SOURCE_DIR}/scripts/output/"
    VERBATIM
    USES_TERMINAL
)

if (USE_OPENMP)
    add_custom_target(run-openmp
        COMMAND sbatch --export=PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR},PROJECT_BINARY_DIR=${CMAKE_RUNTIME_OUTPUT_DIRECTORY},COMPILER_ID=${CMAKE_C_COMPILER_ID} ${PROJECT_SOURCE_DIR}/scripts/jobscript_run.openmp.slurm
        DEPENDS fwi
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/scripts/output
        COMMENT "outputs will be in ${PROJECT_SOURCE_DIR}/scripts/output/"
        VERBATIM
        USES_TERMINAL
    )
endif (USE_OPENMP)

if (USE_OPENACC)
    add_custom_target(run-openacc
        COMMAND sbatch --export=PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR},PROJECT_BINARY_DIR=${CMAKE_RUNTIME_OUTPUT_DIRECTORY},COMPILER_ID=${CMAKE_C_COMPILER_ID} ${PROJECT_SOURCE_DIR}/scripts/jobscript_run.openacc.slurm
        DEPENDS fwi
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/scripts/output
        COMMENT "outputs will be in ${PROJECT_SOURCE_DIR}/scripts/output/"
        VERBATIM
        USES_TERMINAL
    )

    add_custom_target(run-mpi+openacc
        COMMAND sbatch --export=PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR},PROJECT_BINARY_DIR=${CMAKE_RUNTIME_OUTPUT_DIRECTORY},COMPILER_ID=${CMAKE_C_COMPILER_ID} ${PROJECT_SOURCE_DIR}/scripts/jobscript_run.mpi+openacc.slurm
        DEPENDS fwi
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/scripts/output
        COMMENT "outputs will be in ${PROJECT_SOURCE_DIR}/scripts/output/"
        VERBATIM
        USES_TERMINAL
    )
endif (USE_OPENACC)
