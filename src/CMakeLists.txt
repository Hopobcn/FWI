add_library(fwi-core STATIC
    fwi_core.c
    fwi_sched.c
    fwi_common.c
    fwi_kernel.c
    fwi_constants.c
    fwi_propagator.c
)

if (USE_MPI)
    set_target_properties(fwi-core PROPERTIES
        COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}"
        LINK_FLAGS    "${MPI_C_LINK_FLAGS}"
    )
endif (USE_MPI)

target_include_directories(fwi-core PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

if (USE_MPI)
    target_link_libraries(fwi-core
        ${MPI_C_LIBRARIES}
    )
endif (USE_MPI)

if (USE_CUDA_KERNELS)
    add_library(fwi-core-cuda SHARED
        fwi_propagator.cu
    )

    set_target_properties(fwi-core-cuda PROPERTIES
        CUDA_STANDARD 11
    )

    target_include_directories(fwi-core-cuda PUBLIC
        ${PROJECT_SOURCE_DIR}/include
    )

    target_link_libraries(fwi-core
        fwi-core-cuda
    )
endif (USE_CUDA_KERNELS)

if (TRACE_CUDA)
    target_link_libraries(fwi-core
        ${CUDA_TRACE_LIBRARY}
    )
endif (TRACE_CUDA)

add_subdirectory(datagen)
add_subdirectory(schedgen)
