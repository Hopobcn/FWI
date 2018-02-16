cmake_minimum_required(VERSION 3.8)
project(FWI LANGUAGES C)

###### PROJECT OPTIONS #######
option(ENABLE_TESTS     "Enable Unit Tests" OFF)
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

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -c99 -Minfo=accel,mp")

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
    set(FWI_PERFORMANCE_FLAGS "")

    if ("${CMAKE_C_COMPILER_ID}" STREQUAL "PGI")

        set(FWI_PERFORMANCE_FLAGS "-O3 -fast -tp=k8 -Minline")

    elseif (("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU") OR
            ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang"))

        if ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
            set(FWI_PERFORMANCE_FLAGS "-O3 -march=x86-64")
        else ()
            set(FWI_PERFORMANCE_FLAGS "-O3")
        endif ()

    elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")

        set(FWI_PERFORMANCE_FLAGS "-O3 -xHost -restrict")

    else()
        message(WARNING "Compiler not supported!")
    endif ()

    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${FWI_PERFORMANCE_FLAGS}")
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
            set(OpenACC_C_FLAGS "${OpenACC_C_FLAGS} -ta=tesla,cuda8.0,cc20,cc30,cc35,cc50,cc60,lineinfo,maxregcount:128,managed")
        else ()
            #FLAGS FOR POWER8+P100  (Marenostrum4)
            set(OpenACC_C_FLAGS "${OpenACC_C_FLAGS} -ta=tesla,cuda8.0,cc60,maxregcount:128,managed")
        endif ()

        if (USE_CUDA_KERNELS)
            enable_language(CUDA)

            set(CMAKE_CUDA_FLAGS "-Xcompiler \"${COMMON_CXX_GNU_FLAGS}\" -gencode arch=compute_35,code=[sm_35,sm_37] -gencode arch=compute_50,code=sm_52 -gencode arch=compute_60,code=[sm_60,sm_61]")

            set(CMAKE_CUDA_FLAGS_DEBUG   "-O0 -g -G")
            set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -DNDEBUG -res-usage -lineinfo")

            add_definitions("-DUSE_CUDA")

            # list of optimized kernels
            add_definitions("-DVCELL_TL")
            add_definitions("-DVCELL_TR")
            add_definitions("-DVCELL_BR")
            add_definitions("-DVCELL_BL")

            add_definitions("-DSCELL_TR")
            add_definitions("-DSCELL_TL")
            add_definitions("-DSCELL_BR")
            add_definitions("-DSCELL_BL")

            message(STATUS "----[CUDA ENABLED]")

        endif ()

        if (TRACE_CUDA)
            add_definitions("-DTRACE_CUDA")

            list(GET CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES 0 CUDA_TOOLKIT_INCLUDE_DIR)

            find_library(CUDA_TRACE_LIBRARY
                NAMES nvToolsExt
                HINTS "${CUDA_TOOLKIT_INCLUDE_DIR}/../lib64"
                PATHS /usr/local/cuda/lib64
            )
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
message(STATUS "CMAKE_CUDA_FLAGS:               ${CMAKE_CUDA_FLAGS}")
message(STATUS "CMAKE_CUDA_FLAGS_DEBUG:         ${CMAKE_CUDA_FLAGS_DEBUG}")
message(STATUS "CMAKE_CUDA_FLAGS_RELEASE:       ${CMAKE_CUDA_FLAGS_RELEASE}")
message(STATUS "----------------------------------------")

add_subdirectory(src)
add_subdirectory(main)

if (ENABLE_TESTS)
    include(CTest)
    enable_testing()

    add_subdirectory(external)
    add_subdirectory(test)
endif (ENABLE_TESTS)
