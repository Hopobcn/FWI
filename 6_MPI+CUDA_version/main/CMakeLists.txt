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
endif (USE_CUDA_KERNELS)
