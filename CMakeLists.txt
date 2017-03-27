cmake_minimum_required(VERSION 3.0)
project(FWI CXX)

###### PROJECT OPTIONS #######
option(ENABLE_TESTS     "Enable Unit Tests" OFF)
option(PERFORM_IO       "Enable IO"         OFF)
option(IO_STATS         "Print IO stats"    OFF)
option(USE_MPI          "Use MPI"           OFF)
option(USE_CUDA         "Use CUDA"          ON)
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
set(DISABLE_CXX_GNU_WARNINGS "-Wno-unused-parameter")
set(COMMON_CXX_GNU_FLAGS "-Wall -Wextra -Wpointer-arith -Wcast-align -Wwrite-strings -Wswitch-default -Wunreachable-code -Winit-self -Wmissing-field-initializers -Wno-unknown-pragmas -Wundef -Wmissing-declarations ${DISABLE_CXX_GNU_WARNINGS}")


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -c99 -Minfo=accel,inline,mp")

elseif ("${}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ${COMMON_CXX_GNU_FLAGS}")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ${COMMON_CXX_GNU_FLAGS}")


elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c99 -prec-div -fp-model precise -ftz -qopt-prefetch-distance=64 ${COMMON_GNU_FLAGS}")
else()
    message(WARNING "Compiler not supported!")
endif ()
###########################################

###### BUILD-SPECIFIC COMPILE FLAGS #######
if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    # set DEBUG FLAGS (seems that ICC/GCC/PGC understant '-g' and '-O0' great!)
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
    add_definitions("-DDEBUG")
else ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    # set RELEASE FLAGS (unfortunately it's compiler dependent...)
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")

        set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fast -Minline")

    elseif (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") OR
            ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"))

        set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")

    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")

        set(CMAKE_CXX_FLAGS_RELEASE "-O3 -xHost -restrict")

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


if (USE_CUDA)

    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

        find_package(CUDA)

        if (CUDA_FOUND)
            set(CUDA_CXX_FLAGS "-x cuda --cuda-gpu-arch=sm_61 --Xcuda-ptxas='-v' ")
            set(CUDA_EXE_LINKER_FLAGS "--cuda-path=${CUDA_TOOLKIT_ROOT_DIR}")

            include_directories(${CUDA_INCLUDE_DIR})

            message(STATUS "--------------------------------------")
            message(STATUS "----[CUDA ENABLED] :                  ")
            message(STATUS "---- CUDA_VERSION:      ${CUDA_VERSION} (${CUDA_VERSION_STRING})")
            message(STATUS "---- CUDA_TOOLKIT_DIR:  ${CUDA_TOOLKIT_ROOT_DIR}")
            message(STATUS "---- CUDA_INCLUDE_DIRS: ${CUDA_INCLUDE_DIRS}")
            message(STATUS "---- CUDA_LIBRARIES:    ${CUDA_LIBRARIES}")
            message(STATUS "--------                      --------")
            message(STATUS "---- CUDA_COMPILER: ${CMAKE_CXX_COMPILER_ID}")
            message(STATUS "---- CUDA_CXX_FLAGS: ${CUDA_CXX_FLAGS}")
            message(STATUS "---- CUDA_EXE_LINKER_FLAGS: ${CUDA_EXE_LINKER_FLAGS}")
            message(STATUS "--------------------------------------")

            add_compile_options(-x cuda --cuda-gpu-arch=sm_61)

            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CUDA_EXE_LINKER_FLAGS}")

            add_definitions("-DUSE_CUDA")

            if (TRACE_CUDA)
               add_definitions("-DTRACE_CUDA")
               set(CUDA_TRACE_LIBRARY "${CUDA_TOOLKIT_ROOT_DIR}/lib64/libnvToolsExt.so")

               message(STATUS "---- TRACE_CUDA: ON")
               message(STATUS "---- CUDA_TRACE_LIBRARY ${CUDA_TRACE_LIBRARY}")
               message(STATUS "--------------------------------------")
           endif (TRACE_CUDA)

        endif (CUDA_FOUND)


    else ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        message(FATAL_ERROR "Compiler (${CMAKE_CXX_COMPILER_ID}) not supported for CUDA compilation in cpp files. Check your compiler!")
    endif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

endif (USE_CUDA)


if (USE_OPENMP)
    find_package(OpenMP REQUIRED QUIET)

    message(STATUS "--------------------------------------")
    message(STATUS "----[OpenMP ENABLED] setting:         ")
    message(STATUS "---- OpemMP_CXX_FLAGS:    ${OpenMP_CXX_FLAGS}")
    message(STATUS "--------------------------------------")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif (USE_OPENMP)


if (PROFILE)

    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")

        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Mprof")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Mprof")

    elseif (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") OR
            ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"))

        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pg")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -g -pg")

    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")

        set(CMAKE_CXX_FLAGS_RELEASE "-g -p")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -g -p")

    else()
        message(WARNING "Compiler not supported!")
    endif ()

    message(STATUS "Profile flags added!")

endif (PROFILE)


message(STATUS "----------------------------------------")
message(STATUS "BUILD TYPE: ${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_CXX_FLAGS:                  ${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_CXX_FLAGS_DEBUG:            ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "CMAKE_CXX_FLAGS_RELEASE:          ${CMAKE_CXX_FLAGS_RELEASE}")
message(STATUS "CMAKE_EXE_LINKER_FLAGS:         ${CMAKE_EXE_LINKER_FLAGS}")
message(STATUS "CMAKE_EXE_LINKER_FLAGS_DEBUG:   ${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
message(STATUS "CMAKE_EXE_LINKER_FLAGS_RELEASE: ${CMAKE_EXE_LINKER_FLAGS_RELEASE}")
message(STATUS "----------------------------------------")

add_subdirectory(src)
add_subdirectory(main)

if (ENABLE_TESTS)
    include(CTest)
    enable_testing()

    add_subdirectory(external)
    add_subdirectory(test)
endif (ENABLE_TESTS)
