add_executable(fwi
    fwi_main.c
)

target_link_libraries(fwi
    fwi-core
    m
)

add_custom_target(run-seq
    COMMAND mnsubmit ${PROJECT_SOURCE_DIR}/scripts/jobscript_${CMAKE_C_COMPILER_ID}_run.sequential.slurm
    DEPENDS fwi
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/scripts/output
    COMMENT "outputs will be in ../scripts/output/"
    VERBATIM
    USES_TERMINAL
)

if (USE_OPENMP)
    add_custom_target(run-openmp
        COMMAND mnsubmit ${PROJECT_SOURCE_DIR}/scripts/jobscript_${CMAKE_C_COMPILER_ID}_run.openmp.slurm
        DEPENDS fwi
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/scripts/output
        COMMENT "outputs will be in ../scripts/output/"
        VERBATIM
        USES_TERMINAL
    )
endif (USE_OPENMP)

if (USE_OPENACC)
    add_custom_target(run-openacc
        COMMAND mnsubmit ${PROJECT_SOURCE_DIR}/scripts/jobscript_${CMAKE_C_COMPILER_ID}_run.openacc.slurm
        DEPENDS fwi
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/scripts/output
        COMMENT "outputs will be in ../scripts/output/"
        VERBATIM
        USES_TERMINAL
    )

    add_custom_target(run-mpi+openacc
        COMMAND mnsubmit ${PROJECT_SOURCE_DIR}/scripts/jobscript_${CMAKE_C_COMPILER_ID}_run.mpi+openacc.slurm
        DEPENDS fwi
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/scripts/output
        COMMENT "outputs will be in ../scripts/output/"
        VERBATIM
        USES_TERMINAL
    )
endif (USE_OPENACC)
