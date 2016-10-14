#include <cuda_runtime.h>
#include <cstdio>

#include "fwi_propagator.cuh"

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////////////////            HACKATHON TODO:                              //
///////////////////////////////////////////////////////////////////////////////
//
// CUDA:
//    - implement CUDA kernels in this FILE
//    - implement CUDA <--> C glue functions (configure & launch kernels)
//
//    - A) functions like IDX,stencil_X etc.. have to be re-implemented HERE
//             OR
//      B) add __host__ __device__ specifiers in C function and modify CMakeList.txt to compile fwi_propagator.c with nvcc and not PGI (problem--> NVCC does not know what OpenACC is..)
//
//     Implement anything you need
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


void compute_component_vcell_TL_cuda ( /* COMPLETE ME */)
{
    /* Configure & launch CUDA kernel */

    CUDA_CHECK(cudaGetLastError());
};

void compute_component_vcell_TR_cuda ( /* COMPLETE ME */ )
{
    /* Configure & launch CUDA kernel */

    CUDA_CHECK(cudaGetLastError());
};

void compute_component_vcell_BR_cuda ( /* COMPLETE ME */ )
{
    /* Configure & launch CUDA kernel */

    CUDA_CHECK(cudaGetLastError());
};

void compute_component_vcell_BL_cuda ( /* COMPLETE ME */ )
{
    /* Configure & launch CUDA kernel */

    CUDA_CHECK(cudaGetLastError());
};

/* ------------------------------------------------------------------------------ */
/*                                                                                */
/*                               CALCULO DE TENSIONES                             */
/*                                                                                */
/* ------------------------------------------------------------------------------ */

void compute_component_scell_TR_cuda ( /* COMPLETE ME */ )
{
    /* Configure & launch CUDA kernel */

    CUDA_CHECK(cudaGetLastError());
};

void compute_component_scell_TL_cuda ( /* COMPLETE ME */ )
{
    /* Configure & launch CUDA kernel */

    CUDA_CHECK(cudaGetLastError());
};

void compute_component_scell_BR_cuda ( /* COMPLETE ME */ )
{
    /* Configure & launch CUDA kernel */

    CUDA_CHECK(cudaGetLastError());
};


void compute_component_scell_BL_cuda ( /* COMPLETE ME */ )
{
    /* Configure & launch CUDA kernel */

    CUDA_CHECK(cudaGetLastError());
};


