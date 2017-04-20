cmake_minimum_required(VERSION 3.0)
project(FWI C CXX)

###### PROJECT OPTIONS #######
option(ENABLE_TESTS     "Enable Unit Tests" ON)
option(PERFORM_IO       "Enable IO"         OFF)
option(IO_STATS         "Print IO stats"    OFF)
option(USE_MPI          "Use MPI"           OFF)
option(USE_OPENACC      "Use OpenACC 2.5"   OFF)
option(USE_OPENMP       "Use OpenMP"        OFF)
option(USE_CUDA_KERNELS "Use CUDA kernels"  OFF)
option(PROFILE          "Add profiling info" OFF)


###### CMAKE WHERE TO STORE BINARY & LIBS ##########
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)


###### Set default build type #############
if (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif (NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)


###### GENERAL COMPILE FLAGS ##############
set(DISABLE_C_GNU_WARNINGS "-Wno-unused-parameter")
set(COMMON_C_GNU_FLAGS "-Wall -Wextra -Wpointer-arith -Wcast-align -Wwrite-strings -Wswitch-default -Wunreachable-code -Winit-self -Wmissing-field-initializers -Wno-unknown-pragmas -Wundef -Wmissing-declarations ${DISABLE_C_GNU_WARNINGS}")
set(COMMON_CXX_GNU_FLAGS "-Wall -Wextra -Wpointer-arith -Wcast-align -Wwrite-strings -Wswitch-default -Wunreachable-code -Winit-self -Wmissing-field-initializers -Wno-unknown-pragmas -Wundef -Wmissing-declarations")


if ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -c99 -Minfo=accel,inline,mp")

elseif (("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU") OR
        ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang"))

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 ${COMMON_C_GNU_FLAGS}")

elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -prec-div -fp-model precise -ftz -qopt-prefetch-distance=64 ${COMMON_GNU_FLAGS}")
else()
    message(WARNING "Compiler not supported!")
endif ()
###########################################

###### BUILD-SPECIFIC COMPILE FLAGS #######
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    # set DEBUG FLAGS (seems that ICC/GCC/PGC understant '-g' and '-O0' great!)
    set(CMAKE_C_FLAGS_DEBUG "-g -O0")
    add_definitions("-DDEBUG")
else ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    # set RELEASE FLAGS (unfortunately it's compiler dependent...)
    if ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")

        set(CMAKE_C_FLAGS_RELEASE "-O3 -fast -Minline")

    elseif (("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU") OR
            ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang"))

        set(CMAKE_C_FLAGS_RELEASE "-O3")

        if ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
            set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native")
        endif ()

    elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")

        set(CMAKE_C_FLAGS_RELEASE "-O3 -xHost -restrict")

    else()
        message(WARNING "Compiler not supported!")
    endif ()

    add_definitions("-DNDEBUG")

endif ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
##########################################





##########################################
############## OPTIONS ###################
##########################################
if (PERFORM_IO)
    if (IO_STATS)
        add_definitions("-DLOG_IO_STATS")
    endif (IO_STATS)
else (PERFORM_IO)
    add_definitions("-DDO_NOT_PERFORM_IO")
endif (PERFORM_IO)


if (USE_MPI)
    find_package(MPI REQUIRED QUIET)

    message(STATUS "--------------------------------------")
    message(STATUS "----[MPI ENABLED] setting:            ")
    message(STATUS "---- MPI_C_COMPILER:      ${MPI_C_COMPILER}")
    message(STATUS "---- MPI_C_COMPILE_FLAGS: ${MPI_C_COMPILE_FLAGS}")
    message(STATUS "---- MPI_C_INCLUDE_PATH:  ${MPI_C_INCLUDE_PATH}")
    message(STATUS "---- MPI_C_LINK_FLAGS:    ${MPI_C_LINK_FLAGS}")
    message(STATUS "---- MPI_C_LIBRARIES:     ${MPI_C_LIBRARIES}")
    message(STATUS "--------------------------------------")

    include_directories(${MPI_C_INCLUDE_PATH})
    add_definitions("-DUSE_MPI")

endif (USE_MPI)


if (USE_OPENACC)

    if ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
        set(OpenACC_C_FLAGS "-acc -ta=multicore")

        if ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
            #Flags for x86 machines (Minotauro, ironman, thor, ...)
            set(OpenACC_C_FLAGS "${OpenACC_C_FLAGS} -ta=tesla,cuda8.0,cc20,cc35,cc50,cc60,lineinfo,managed")
        else ()
            #FLAGS FOR POWER8+P100  (Marenostrum4)
            set(OpenACC_C_FLAGS "${OpenACC_C_FLAGS} -ta=tesla,cuda8.0,cc60,managed")
        endif ()

        find_package(CUDA)

        if (USE_CUDA_KERNELS AND CUDA_FOUND)
            # gcc version should be > 4.8 and < 5.0 (CUDA 8)
            if ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
                set(CUDA_HOST_COMPILER "gcc-5")
            else ()
                set(CUDA_HOST_COMPILER "/usr/bin/gcc")
            endif ()

           # set(CUDA_HOST_COMPILATION_CPP OFF)
            set(CUDA_PROPAGATE_HOST_FLAGS OFF)
            set(CUDA_SEPARABLE_COMPILATION OFF)

            set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11 -Xcompiler \"${COMMON_CXX_GNU_FLAGS}\" -gencode arch=compute_30,code=\"sm_35,sm_37\" -gencode arch=compute_50,code=sm_52 -gencode arch=compute_30,code=\"sm_60,sm_61\" -lineinfo")
            set(CUDA_NVCC_FLAGS_DEBUG   "-O0 -g -G")
            set(CUDA_NVCC_FLAGS_RELEASE "-O3       -res-usage")

            add_definitions("-DUSE_CUDA")
            include_directories(${CUDA_INCLUDE_DIR})

            # list of optimized kernels
            add_definitions("-DVCELL_TL")
            add_definitions("-DVCELL_TR")
            add_definitions("-DVCELL_BR")
            add_definitions("-DVCELL_BL")

            add_definitions("-DSCELL_TR")
            add_definitions("-DSCELL_TL")
            add_definitions("-DSCELL_BR")
            add_definitions("-DSCELL_BL")

            message(STATUS "--------------------------------------")
            message(STATUS "----[CUDA ENABLED] :                  ")
            message(STATUS "---- CUDA_VERSION:      ${CUDA_VERSION} (${CUDA_VERSION_STRING})")
            message(STATUS "---- CUDA_TOOLKIT_DIR:  ${CUDA_TOOLKIT_ROOT_DIR}")
            message(STATUS "---- CUDA_INCLUDE_DIRS: ${CUDA_INCLUDE_DIRS}")
            message(STATUS "---- CUDA_LIBRARIES:    ${CUDA_LIBRARIES}")
            message(STATUS "--------                      --------")
            message(STATUS "---- CUDA_HOST_COMPILER: ${CUDA_HOST_COMPILER}")
            message(STATUS "--------------------------------------")

        endif (USE_CUDA_KERNELS AND CUDA_FOUND)

        if (TRACE_CUDA)
            add_definitions("-DTRACE_CUDA")
            set(CUDA_TRACE_LIBRARY "${CUDA_TOOLKIT_ROOT_DIR}/lib64/libnvToolsExt.so")
        endif (TRACE_CUDA)

    elseif( ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU") AND
            ( ${CMAKE_C_COMPILER_VERSION} VERSION_GREATER 6.1))

        set(OpenACC_C_FLAGS "-fopenacc")


        message(WARNING "GCC has very EXPERIMENTAL support for OpenACC 2.0a. We recommend using PGI 16.5 for OpenACC!")
    else ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")
        message(FATAL_ERROR "Compiler (${CMAKE_C_COMPILER_ID}, ${CMAKE_C_COMPILER_VERSION}) not supported for OpenACC compilation. Check your compiler!")
    endif ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")

    message(STATUS "--------------------------------------")
    message(STATUS "----[OpenACC ENABLED] setting:        ")
    message(STATUS "---- OpenACC_C_FLAGS:   ${OpenACC_C_FLAGS}")
    message(STATUS "--------------------------------------")

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenACC_C_FLAGS}")
endif (USE_OPENACC)


if (USE_OPENMP AND NOT USE_OPENACC)
    find_package(OpenMP REQUIRED QUIET)

    message(STATUS "--------------------------------------")
    message(STATUS "----[OpenMP ENABLED] setting:         ")
    message(STATUS "---- OpemMP_C_FLAGS:    ${OpenMP_C_FLAGS}")
    message(STATUS "--------------------------------------")

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
else (USE_OPENMP AND NOT USE_OPENACC)
    if (USE_OPENMP)
        message(WARNING "OpenMP disabled since OpenACC is enabled!. Disable OpenACC if you want to compile with OpenMP support.")
        set(USE_OPENMP FALSE)
    endif (USE_OPENMP)
endif (USE_OPENMP AND NOT USE_OPENACC)


if (PROFILE)

    if ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")

        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Mprof")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Mprof")

    elseif (("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU") OR
            ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang"))

        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -pg")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -g -pg")

    elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")

        set(CMAKE_C_FLAGS_RELEASE "-g -p")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -g -p")

    else()
        message(WARNING "Compiler not supported!")
    endif ()

    message(STATUS "Profile flags added!")

endif (PROFILE)


message(STATUS "----------------------------------------")
message(STATUS "CPU TYPE: ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "BUILD TYPE: ${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_C_FLAGS:                  ${CMAKE_C_FLAGS}")
message(STATUS "CMAKE_C_FLAGS_DEBUG:            ${CMAKE_C_FLAGS_DEBUG}")
message(STATUS "CMAKE_C_FLAGS_RELEASE:          ${CMAKE_C_FLAGS_RELEASE}")
message(STATUS "CMAKE_EXE_LINKER_FLAGS:         ${CMAKE_EXE_LINKER_FLAGS}")
message(STATUS "CMAKE_EXE_LINKER_FLAGS_DEBUG:   ${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
message(STATUS "CMAKE_EXE_LINKER_FLAGS_RELEASE: ${CMAKE_EXE_LINKER_FLAGS_RELEASE}")
message(STATUS "CUDA_NVCC_FLAGS:                ${CUDA_NVCC_FLAGS}")
message(STATUS "CUDA_NVCC_FLAGS_DEBUG:          ${CUDA_NVCC_FLAGS_DEBUG}")
message(STATUS "CUDA_NVCC_FLAGS_RELEASE:        ${CUDA_NVCC_FLAGS_RELEASE}")
message(STATUS "----------------------------------------")

add_subdirectory(src)
add_subdirectory(main)

if (ENABLE_TESTS)
    include(CTest)
    enable_testing()

    add_subdirectory(external)
    add_subdirectory(test)
endif (ENABLE_TESTS)
