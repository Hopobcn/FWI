add_library(fwi-core STATIC
    fwi_core.cpp
    fwi_common.cpp
    fwi_kernel.cpp
    fwi_propagator.cpp
)

if (USE_MPI)
    set_target_properties(fwi-core PROPERTIES
        COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}"
        LINK_FLAGS    "${MPI_CXX_LINK_FLAGS}"
    )
endif (USE_MPI)

target_include_directories(fwi-core PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

if (USE_MPI)
    target_link_libraries(fwi-core
        ${MPI_CXX_LIBRARIES}
    )
endif (USE_MPI)

#f (USE_CUDA_KERNELS)
#   cuda_include_directories( ${PROJECT_SOURCE_DIR}/include )
#
#   cuda_add_library(fwi-core-cuda SHARED
#       fwi_common.cu
#       fwi_propagator.cu
#   )
#
#   target_link_libraries(fwi-core-cuda
#       ${CUDA_LIBRARIES}
#   )
#
#   target_link_libraries(fwi-core
#       fwi-core-cuda
#   )
#ndif (USE_CUDA_KERNELS)
#
#f (TRACE_CUDA)
#   target_link_libraries(fwi-core
#       ${CUDA_TRACE_LIBRARY}
#   )
#ndif (TRACE_CUDA)

add_subdirectory(datagen)
