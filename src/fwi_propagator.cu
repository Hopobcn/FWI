#include <cuda_runtime.h>
#include <cstdio>

#include "fwi_propagator.cuh"

//// IMPLEMENT ANYTHING YOU NEED
// You could either:
//  - try to modify host functions present in fwi_propagator.c to be __host__ __device__ functions
//     (you will have to deal with C<-->C++ interoperability (.cu files are treated as C++)
//  - or implement your own specialization here


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


