cmake_minimum_required(VERSION 3.0)
project(FWI C)

option(ENABLE_TESTS     "Enable Unit Tests" ON)
option(PERFORM_IO       "Enable IO"         OFF)
option(IO_STATS         "Print IO stats"    OFF)
option(USE_MPI          "Use MPI"           OFF)
option(USE_OPENMP       "Use OpenMP"        OFF)
option(USE_OPENACC      "Use OpenACC 2.5"   OFF)
option(USE_CUDA_KERNELS "Use CUDA kernels"  OFF)


###### GENERAL COMPILE FLAGS ##############
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
    
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -c99")

elseif (("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU") OR
        ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang"))

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -Wall")

elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -prec-div -fp-model precise -ftz -qopt-prefetch-distance=64")
else()
    message(WARNING "Compiler not supported!")
endif ()
###########################################

###### BUILD-SPECIFIC COMPILE FLAGS #######
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    # set DEBUG FLAGS (seems that ICC/GCC/PGC understant '-g' and '-O0' great!)
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -O0")
    add_definitions("-DDEBUG")
else ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    # set RELEASE FLAGS (unfortunately it's compiler dependent...)
    if ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
       
       set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fast")
    
    elseif (("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU") OR
           ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang"))
    
       set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native")
    
    elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")

        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -xHost -restrict")

    else()
        message(WARNING "Compiler not supported!")
    endif ()

endif ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
##########################################





##########################################
############## OPTIONS ###################
##########################################
if (PERFORM_IO)
    if (IO_STATS)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DLOG_IO_STATS")
    endif (IO_STATS)
else (PERFORM_IO)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDO_NOT_PERFORM_IO")
endif (PERFORM_IO)

if (USE_MPI)
    find_package(MPI REQUIRED QUIET)

    message(STATUS "MPI_C_COMPILER:       ${MPI_C_COMPILER}")
    message(STATUS "MPI_C_COMPILE_FLAGS:  ${MPI_C_COMPILE_FLAGS}")
    message(STATUS "MPI_C_INCLUDE_PATH:   ${MPI_C_INCLUDE_PATH}")
    message(STATUS "MPI_C_LINK_FLAGS:     ${MPI_C_LINK_FLAGS}")
    message(STATUS "MPI_C_LIBRARIES:      ${MPI_C_LIBRARIES}")

    include_directories(${MPI_INCLUDE_PATH})
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS}")
    #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_C_LINK_FLAGS}")
endif (USE_MPI)

if (USE_OPENMP)
    find_package(OpenMP REQUIRED QUIET)

    message(STATUS "OpemMP_C_FLAGS ${OpenMP_C_FLAGS}")

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif (USE_OPENMP)

if (USE_OPENACC)

    if ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -acc -ta=tesla,cc20,cc35,cc50,maxregcount:128,lineinfo")

        if (USE_CUDA_KERNELS)
            find_package(CUDA REQUIRED)

            set(CUDA_HOST_COMPILER "gcc")

            add_definitions("-DUSE_CUDA")
            include_directories(${CUDA_INCLUDE_DIR})
        endif (USE_CUDA_KERNELS)
    elseif( ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU") AND
            ( ${CMAKE_C_COMPILER_VERSION} VERSION_GREATER 6.1))
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenacc")

        message(WARNING "GCC has very EXPERIMENTAL support for OpenACC 2.0a. We recommend using PGI 16.5 for OpenACC!")
    else ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
        message(FATAL_ERROR "Compiler (${CMAKE_C_COMPILER_ID}, ${CMAKE_C_COMPILER_VERSION}) not supported for OpenACC compilation. Check your compiler!")
    endif ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")

endif (USE_OPENACC)


add_subdirectory(src)

if (ENABLE_TESTS)
    include(CTest)
    enable_testing()

    add_subdirectory(external)
    add_subdirectory(test)
endif (ENABLE_TESTS)
