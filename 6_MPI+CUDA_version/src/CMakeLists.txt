add_library(fwi-core STATIC
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
    cuda_add_library(fwi-core-cuda STATIC
        fwi_common.cu
        fwi_propagator.cu
    )
endif (USE_CUDA_KERNELS)

add_executable(fwi
    fwi_main.c
)

target_link_libraries(fwi
    fwi-core
    m
)

if (USE_CUDA_KERNELS)
    target_link_libraries(fwi
        fwi-core-cuda
    )
endif ()
