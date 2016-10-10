add_library(fwi-core STATIC
    fwi_core.c
    fwi_common.c
    fwi_kernel.c
    fwi_propagator.c
)

if (USE_MPI)
    set_target_properties(fwi-core PROPERTIES
        COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}"
        LINK_FLAGS    "${MPI_C_LINK_FLAGS}"
    )
endif (USE_MPI)

target_include_directories(fwi-core PUBLIC
    ${PROJECT_SOURCE_DIR}/include/fwi
)

if (USE_MPI)
    target_link_libraries(fwi-core
        ${MPI_C_LIBRARIES}
    )
endif (USE_MPI)

if (USE_CUDA_KERNELS)
    cuda_include_directories( ${PROJECT_SOURCE_DIR}/include/fwi )

    cuda_add_library(fwi-core-cuda SHARED
        fwi_propagator.cu
    )

    target_link_libraries(fwi-core-cuda
        ${CUDA_LIBRARIES}
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
