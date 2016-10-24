/*
 * =============================================================================
 * Copyright (c) 2016, Barcelona Supercomputing Center (BSC)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * =============================================================================
 */

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


