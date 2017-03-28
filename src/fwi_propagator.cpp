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

#include "fwi/fwi_propagator.hpp"

HOST_DEVICE_INLINE
integer IDX (const integer z,
             const integer x,
             const integer y,
             const integer dimmz,
             const integer dimmx)
{
    return ((y*dimmx)+x)*dimmz + z;
};

HOST_DEVICE_INLINE
real stencil_Z ( const integer off,
                 const real*   ptr,
                 const real    dzi,
                 const integer z,
                 const integer x,
                 const integer y,
                 const integer dimmz,
                 const integer dimmx)
{
    return  ((C0 * ( ptr[IDX(z  +off,x,y,dimmz,dimmx)] - ptr[IDX(z-1+off,x,y,dimmz,dimmx)]) +
              C1 * ( ptr[IDX(z+1+off,x,y,dimmz,dimmx)] - ptr[IDX(z-2+off,x,y,dimmz,dimmx)]) +
              C2 * ( ptr[IDX(z+2+off,x,y,dimmz,dimmx)] - ptr[IDX(z-3+off,x,y,dimmz,dimmx)]) +
              C3 * ( ptr[IDX(z+3+off,x,y,dimmz,dimmx)] - ptr[IDX(z-4+off,x,y,dimmz,dimmx)])) * dzi );
};

template <const int HALO,
          const int BDIMX>
DEVICE inline
real stencil_Z_shfl (const integer off,
                     const real* __restrict__ ptr_gmem,
                     const real    dzi,
                     const integer z,
                     const integer x,
                     const integer y,
                     const integer dimmz,
                     const integer dimmx)
{
    real current = (z+off < dimmz && x < dimmx) ? ptr_gmem[IDX(z+off,x,y,dimmz,dimmx)] : 0.0f;
    real right3  = __shfl_down(current, 3);
    real right2  = __shfl_down(current, 2);
    real right1  = __shfl_down(current, 1);
    real left1   = __shfl_up(current, 1);
    real left2   = __shfl_up(current, 2);
    real left3   = __shfl_up(current, 3);
    real left4   = __shfl_up(current, 4);

    /* For threads without neighbors: */
    if (threadIdx.x < 1 /* 1 */) left1 = (z+off-1 < dimmz && x < dimmx) ? ptr_gmem[IDX(z+off-1,x,y,dimmz,dimmx)] : 0.0f;
    if (threadIdx.x < 2 /* 2 */) left2 = (z+off-2 < dimmz && x < dimmx) ? ptr_gmem[IDX(z+off-2,x,y,dimmz,dimmx)] : 0.0f;
    if (threadIdx.x < 3 /* 3 */) left3 = (z+off-3 < dimmz && x < dimmx) ? ptr_gmem[IDX(z+off-3,x,y,dimmz,dimmx)] : 0.0f;
    if (threadIdx.x < 4 /* 4 */) left4 = (z+off-4 < dimmz && x < dimmx) ? ptr_gmem[IDX(z+off-4,x,y,dimmz,dimmx)] : 0.0f;

    if (threadIdx.x >= BDIMX-1 /* 1 */) right1 = (z+off+1 < dimmz && x < dimmx) ? ptr_gmem[IDX(z+off+1,x,y,dimmz,dimmx)] : 0.0f;
    if (threadIdx.x >= BDIMX-2 /* 2 */) right2 = (z+off+1 < dimmz && x < dimmx) ? ptr_gmem[IDX(z+off+2,x,y,dimmz,dimmx)] : 0.0f;
    if (threadIdx.x >= BDIMX-3 /* 3 */) right3 = (z+off+1 < dimmz && x < dimmx) ? ptr_gmem[IDX(z+off+3,x,y,dimmz,dimmx)] : 0.0f;

    return  ((C0 * ( current - left1) +
              C1 * ( right1  - left2) +
              C2 * ( right2  - left3) +
              C3 * ( right3  - left4)) * dzi );
};

HOST_DEVICE_INLINE
real stencil_X( const integer off,
                const real*   ptr,
                const real    dxi,
                const integer z,
                const integer x,
                const integer y,
                const integer dimmz,
                const integer dimmx)
{
    return ((C0 * ( ptr[IDX(z,x  +off,y,dimmz,dimmx)] - ptr[IDX(z,x-1+off,y,dimmz,dimmx)]) +
             C1 * ( ptr[IDX(z,x+1+off,y,dimmz,dimmx)] - ptr[IDX(z,x-2+off,y,dimmz,dimmx)]) +
             C2 * ( ptr[IDX(z,x+2+off,y,dimmz,dimmx)] - ptr[IDX(z,x-3+off,y,dimmz,dimmx)]) +
             C3 * ( ptr[IDX(z,x+3+off,y,dimmz,dimmx)] - ptr[IDX(z,x-4+off,y,dimmz,dimmx)])) * dxi );
};

template <const int NX,
          const int NY,
          const int HALO,
          const int BDIMY>
DEVICE inline
real stencil_X_smem(const integer off,
                          real    ptr_smem[NY][NX],
                    const real* __restrict__ ptr_gmem,
                    const real    dxi,
                    const integer z,
                    const integer x,
                    const integer y,
                    const integer dimmz,
                    const integer dimmx)
{
    //__syncthreads();
    const integer tx = threadIdx.x;
    const integer ty = threadIdx.y+HALO;
    ///////////// intra-block communication///////////////////
    ptr_smem[ty][tx] = ((z < dimmz && x+off < dimmx) ? ptr_gmem[IDX(z,x+off,y,dimmz,dimmx)] : 0.0f);
    if (threadIdx.y < HALO)
    {
        ptr_smem[ty-HALO ][tx] = ((z < dimmz && x+off-HALO  < dimmx) ? ptr_gmem[IDX(z,x+off-HALO, y,dimmz,dimmx)] : 0.0f);
        ptr_smem[ty+BDIMY][tx] = ((z < dimmz && x+off+BDIMY < dimmx) ? ptr_gmem[IDX(z,x+off+BDIMY,y,dimmz,dimmx)] : 0.0f);
    }
    __syncthreads();

    return ((C0 * ( ptr_smem[ty  ][tx] - ptr_smem[ty-1][tx] ) +
             C1 * ( ptr_smem[ty+1][tx] - ptr_smem[ty-2][tx] ) +
             C2 * ( ptr_smem[ty+2][tx] - ptr_smem[ty-3][tx] ) +
             C3 * ( ptr_smem[ty+3][tx] - ptr_smem[ty-4][tx] )) * dxi );
};

HOST_DEVICE_INLINE
real stencil_Y( const integer off,
                const real*   ptr,
                const real    dyi,
                const integer z,
                const integer x,
                const integer y,
                const integer dimmz,
                const integer dimmx)
{
    return ((C0 * ( ptr[IDX(z,x,y  +off,dimmz,dimmx)] - ptr[IDX(z,x,y-1+off,dimmz,dimmx)]) +
             C1 * ( ptr[IDX(z,x,y+1+off,dimmz,dimmx)] - ptr[IDX(z,x,y-2+off,dimmz,dimmx)]) +
             C2 * ( ptr[IDX(z,x,y+2+off,dimmz,dimmx)] - ptr[IDX(z,x,y-3+off,dimmz,dimmx)]) +
             C3 * ( ptr[IDX(z,x,y+3+off,dimmz,dimmx)] - ptr[IDX(z,x,y-4+off,dimmz,dimmx)])) * dyi );
};

/* -------------------------------------------------------------------- */
/*                     KERNELS FOR VELOCITY                             */
/* -------------------------------------------------------------------- */

HOST_DEVICE_INLINE
real rho_BL ( const real*   rho,
              const integer z,
              const integer x,
              const integer y,
              const integer dimmz,
              const integer dimmx)
{
    return (2.0f / (rho[IDX(z,x,y,dimmz,dimmx)] + rho[IDX(z+1,x,y,dimmz,dimmx)]));
};

template <const int BDIMX>
DEVICE inline
real rho_BL_shfl (const real* __restrict__ ptr_gmem,
                  const integer z,
                  const integer x,
                  const integer y,
                  const integer dimmz,
                  const integer dimmx)
{
    real current = (z<dimmz && x<dimmx) ? ptr_gmem[IDX(z,x,y,dimmz,dimmx)] : 0.0f;
    real right1  = __shfl_down(current, 1);

    /* For threads without neighbors: */
    if (threadIdx.x >= BDIMX-1 /* 1 */) right1 = (z+1<dimmz && x<dimmx) ? ptr_gmem[IDX(z+1,x,y,dimmz,dimmx)] : 0.0f;

    return (2.0f / (current + right1));
};

HOST_DEVICE_INLINE
real rho_TR ( const real*   rho,
              const integer z,
              const integer x,
              const integer y,
              const integer dimmz,
              const integer dimmx)
{
    return (2.0f / (rho[IDX(z,x,y,dimmz,dimmx)] + rho[IDX(z,x+1,y,dimmz,dimmx)]));
};

template <const int BDIMX,
          const int BDIMY>
DEVICE inline
real rho_TR_smem (      real rho_smem[BDIMY+1][BDIMX],
                  const real* __restrict__ rho_gmem,
                  const integer z,
                  const integer x,
                  const integer y,
                  const integer dimmz,
                  const integer dimmx)
{
    const integer tx = threadIdx.x;
    const integer ty = threadIdx.y;

    rho_smem[ty][tx] = (z<dimmz && x<dimmx) ? rho_gmem[IDX(z,x,y,dimmz,dimmx)] : 0.0f;

    /* For threads without neighbors: */
    if (ty < 1)
        rho_smem[ty+BDIMY][tx] = (z<dimmz && x+BDIMY<dimmx) ? rho_gmem[IDX(z,x+BDIMY,y,dimmz,dimmx)] : 0.0f;

    __syncthreads();

    return (2.0f / (rho_smem[ty][tx] + rho_smem[ty+1][tx]));
};

HOST_DEVICE_INLINE
real rho_BR ( const real*   rho,
              const integer z,
              const integer x,
              const integer y,
              const integer dimmz,
              const integer dimmx)
{
    return ( 8.0f/ ( rho[IDX(z  ,x  ,y  ,dimmz,dimmx)] +
                     rho[IDX(z+1,x  ,y  ,dimmz,dimmx)] +
                     rho[IDX(z  ,x+1,y  ,dimmz,dimmx)] +
                     rho[IDX(z  ,x  ,y+1,dimmz,dimmx)] +
                     rho[IDX(z  ,x+1,y+1,dimmz,dimmx)] +
                     rho[IDX(z+1,x+1,y  ,dimmz,dimmx)] +
                     rho[IDX(z+1,x  ,y+1,dimmz,dimmx)] +
                     rho[IDX(z+1,x+1,y+1,dimmz,dimmx)]) );
};

template <const int BDIMX,
          const int BDIMY>
DEVICE inline
real rho_BR_smem (      real rho_smem[BDIMY+1][BDIMX+1],
                  const real* __restrict__ rho_gmem,
                  const real rho_current,
                  const real rho_front1,
                  const integer z,
                  const integer x,
                  const integer y,
                  const integer dimmz,
                  const integer dimmx)
{
    const integer tx = threadIdx.x;
    const integer ty = threadIdx.y;

    rho_smem[ty][tx] = (z<dimmz && x<dimmx) ? rho_gmem[IDX(z,x,y,dimmz,dimmx)] : 0.0f;

    /* For threads without neighbors: */
    if (ty < 1)
        rho_smem[ty+BDIMY][tx      ] = (z<dimmz && x+BDIMY<dimmx) ? rho_gmem[IDX(z,x+BDIMY,y,dimmz,dimmx)] : 0.0f;
    if (tx < 1)
        rho_smem[ty      ][tx+BDIMX] = (z+BDIMX<dimmz && x<dimmx) ? rho_gmem[IDX(z+BDIMX,x,y,dimmz,dimmx)] : 0.0f;

    __syncthreads();

    return (8.0f/ ( rho_current                            +
                    rho_smem[ty  ][tx+1]                   +
                    rho_smem[ty+1][tx  ]                   +
                    rho_front1                             +
                    rho_gmem[IDX(z  ,x+1,y+1,dimmz,dimmx)] +
                    rho_gmem[IDX(z+1,x+1,y  ,dimmz,dimmx)] +
                    rho_gmem[IDX(z+1,x  ,y+1,dimmz,dimmx)] +
                    rho_gmem[IDX(z+1,x+1,y+1,dimmz,dimmx)]) );
};

HOST_DEVICE_INLINE
real rho_TL ( const real*   rho,
              const integer z,
              const integer x,
              const integer y,
              const integer dimmz,
              const integer dimmx)
{
    return (2.0f / (rho[IDX(z,x,y,dimmz,dimmx)] + rho[IDX(z,x,y+1,dimmz,dimmx)]));
};

//////////////////////////////////// VCELL TL /////////////////////////////////////

template <const int HALO  = 4,
          const int BDIMX = 32,
          const int BDIMY = 4,
          const int NX    = BDIMX+2*HALO,
          const int NY    = BDIMY+2*HALO>
__global__
__launch_bounds__(128, 8)
void compute_component_vcell_TL_cuda ( real* __restrict__ vptr,
                                 const real* __restrict__ szptr,
                                 const real* __restrict__ sxptr,
                                 const real* __restrict__ syptr,
                                 const real* __restrict__ rho,
                                 const real               dt,
                                 const real               dzi,
                                 const real               dxi,
                                 const real               dyi,
                                 const integer            nz0,
                                 const integer            nzf,
                                 const integer            nx0,
                                 const integer            nxf,
                                 const integer            ny0,
                                 const integer            nyf,
                                 const offset_t           SZ,
                                 const offset_t           SX,
                                 const offset_t           SY,
                                 const integer            dimmz,
                                 const integer            dimmx)
{
    integer z = blockIdx.x * blockDim.x + threadIdx.x + nz0;
    integer x = blockIdx.y * blockDim.y + threadIdx.y + nx0;

    // WARNING: We can't predicate threads that fall outside of the [nz0:nzf][nx0:nxf] range because
    //          we use COLLECTIVE operations like SHUFFLE & SHARED.
    //          PREVENT incorrect GMEM memory access by checking boundaries at every access

    __shared__ real sx_smem[NY][NX];
    real sy_front1, sy_front2, sy_front3;
    real sy_back1, sy_back2, sy_back3, sy_back4;
    real sy_current;

    sy_back3   = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+0+SY,dimmz,dimmx)] : 0.0f;
    sy_back2   = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+1+SY,dimmz,dimmx)] : 0.0f;
    sy_back1   = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+2+SY,dimmz,dimmx)] : 0.0f;
    sy_current = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+3+SY,dimmz,dimmx)] : 0.0f;
    sy_front1  = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+4+SY,dimmz,dimmx)] : 0.0f;
    sy_front2  = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+5+SY,dimmz,dimmx)] : 0.0f;
    sy_front3  = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+6+SY,dimmz,dimmx)] : 0.0f;

    real rho_current, rho_front1;
    rho_front1 = (z<dimmz && x<dimmx) ? rho[IDX(z,x,ny0,dimmz,dimmx)] : 0.0f;

    for(integer y = ny0;
                y < nyf;
                y++)
    {
        /////////// register tiling-advance plane ////////////////
        sy_back4   = sy_back3;
        sy_back3   = sy_back2;
        sy_back2   = sy_back1;
        sy_back1   = sy_current;
        sy_current = sy_front1;
        sy_front1  = sy_front2;
        sy_front2  = sy_front3;
        sy_front3  = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,y+SY+HALO-1,dimmz,dimmx)] : 0.0f;
        ///////////////////////
        rho_current = rho_front1;
        rho_front1  = (z<dimmz && x<dimmx) ? rho[IDX(z,x,y+1,dimmz,dimmx)] : 0.0f;
        //////////////////////////////////////////////////////////

        const real lrho = (2.0f / (rho_current + rho_front1));

        const real stz = stencil_Z_shfl <HALO,BDIMX> (SZ, szptr, dzi, z, x, y, dimmz, dimmx);

        const real stx = stencil_X_smem <NX,NY,HALO,BDIMY> (SX, sx_smem, sxptr, dxi, z, x, y, dimmz, dimmx);

        const real sty = ((C0 * ( sy_current - sy_back1 ) +
                            C1 * ( sy_front1  - sy_back2 ) +
                            C2 * ( sy_front2  - sy_back3 ) +
                            C3 * ( sy_front3  - sy_back4 )) * dyi );

        if (z < nzf && x < nxf)
        {
            vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
        }
    }
}

__global__
void compute_component_vcell_TL_cuda ( real*    vptr,
                                 const real*    szptr,
                                 const real*    sxptr,
                                 const real*    syptr,
                                 const real*    rho,
                                 const real     dt,
                                 const real     dzi,
                                 const real     dxi,
                                 const real     dyi,
                                 const integer  nz0,
                                 const integer  nzf,
                                 const integer  nx0,
                                 const integer  nxf,
                                 const integer  ny0,
                                 const integer  nyf,
                                 const offset_t SZ,
                                 const offset_t SX,
                                 const offset_t SY,
                                 const integer  dimmz,
                                 const integer  dimmx)
{
    for(integer z = blockIdx.x * blockDim.x + threadIdx.x + nz0;
                z < nzf;
                z += gridDim.x * blockDim.x)
    {
        for(integer x = blockIdx.y * blockDim.y + threadIdx.y + nx0;
                    x < nxf;
                    x += gridDim.y * blockDim.y)
        {
            for(integer y = ny0;
                        y < nyf;
                        y++)
            {
                const real lrho = rho_TL(rho, z, x, y, dimmz, dimmx);

                const real stx  = stencil_X( SX, sxptr, dxi, z, x, y, dimmz, dimmx);
                const real sty  = stencil_Y( SY, syptr, dyi, z, x, y, dimmz, dimmx);
                const real stz  = stencil_Z( SZ, szptr, dzi, z, x, y, dimmz, dimmx);

                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
}

void compute_component_vcell_TL (      real*    vptr,
                                 const real*    szptr,
                                 const real*    sxptr,
                                 const real*    syptr,
                                 const real*    rho,
                                 const real     dt,
                                 const real     dzi,
                                 const real     dxi,
                                 const real     dyi,
                                 const integer  nz0,
                                 const integer  nzf,
                                 const integer  nx0,
                                 const integer  nxf,
                                 const integer  ny0,
                                 const integer  nyf,
                                 const offset_t SZ,
                                 const offset_t SX,
                                 const offset_t SY,
                                 const integer  dimmz,
                                 const integer  dimmx)
{
#if defined(_OPENMP)
    #pragma omp parallel for
#endif /* end _OPENACC */
    for(integer y=ny0; y < nyf; y++)
    {
        for(integer x=nx0; x < nxf; x++)
        {
            for(integer z=nz0; z < nzf; z++)
            {
                const real lrho = rho_TL(rho, z, x, y, dimmz, dimmx);

                const real stx  = stencil_X(SX, sxptr, dxi, z, x, y, dimmz, dimmx);
                const real sty  = stencil_Y(SY, syptr, dyi, z, x, y, dimmz, dimmx);
                const real stz  = stencil_Z(SZ, szptr, dzi, z, x, y, dimmz, dimmx);

                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
};


//////////////////////////////////// VCELL TR /////////////////////////////////////

template <const int HALO  = 4,
          const int BDIMX = 32,
          const int BDIMY = 4,
          const int NX    = BDIMX+2*HALO,
          const int NY    = BDIMY+2*HALO>
__global__
__launch_bounds__(128, 8)
void compute_component_vcell_TR_cuda ( real* __restrict__ vptr,
                                 const real* __restrict__ szptr,
                                 const real* __restrict__ sxptr,
                                 const real* __restrict__ syptr,
                                 const real* __restrict__ rho,
                                 const real               dt,
                                 const real               dzi,
                                 const real               dxi,
                                 const real               dyi,
                                 const integer            nz0,
                                 const integer            nzf,
                                 const integer            nx0,
                                 const integer            nxf,
                                 const integer            ny0,
                                 const integer            nyf,
                                 const offset_t           SZ,
                                 const offset_t           SX,
                                 const offset_t           SY,
                                 const integer            dimmz,
                                 const integer            dimmx)
{
    integer z = blockIdx.x * blockDim.x + threadIdx.x + nz0;
    integer x = blockIdx.y * blockDim.y + threadIdx.y + nx0;

    // WARNING: We can't predicate threads that fall outside of the [nz0:nzf][nx0:nxf] range because
    //          we use COLLECTIVE operations like SHUFFLE & SHARED.
    //          PREVENT incorrect GMEM memory access by checking boundaries at every access

    __shared__ real sx_smem[NY][NX];
    real sy_front1, sy_front2, sy_front3;
    real sy_back1, sy_back2, sy_back3, sy_back4;
    real sy_current;

    sy_back3   = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+0+SY,dimmz,dimmx)] : 0.0f;
    sy_back2   = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+1+SY,dimmz,dimmx)] : 0.0f;
    sy_back1   = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+2+SY,dimmz,dimmx)] : 0.0f;
    sy_current = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+3+SY,dimmz,dimmx)] : 0.0f;
    sy_front1  = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+4+SY,dimmz,dimmx)] : 0.0f;
    sy_front2  = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+5+SY,dimmz,dimmx)] : 0.0f;
    sy_front3  = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+6+SY,dimmz,dimmx)] : 0.0f;

    __shared__ real rho_smem[BDIMY+1][BDIMX];

    for(integer y = ny0;
                y < nyf;
                y++)
    {
        /////////// register tiling-advance plane ////////////////
        sy_back4   = sy_back3;
        sy_back3   = sy_back2;
        sy_back2   = sy_back1;
        sy_back1   = sy_current;
        sy_current = sy_front1;
        sy_front1  = sy_front2;
        sy_front2  = sy_front3;
        sy_front3  = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,y+SY+HALO-1,dimmz,dimmx)] : 0.0f;
        //////////////////////////////////////////////////////////

        const real lrho = rho_TR_smem <BDIMX,BDIMY> (rho_smem, rho, z, x, y, dimmz, dimmx);

        const real stz = stencil_Z_shfl <HALO,BDIMX> (SZ, szptr, dzi, z, x, y, dimmz, dimmx);

        const real stx = stencil_X_smem <NX,NY,HALO,BDIMY> (SX, sx_smem, sxptr, dxi, z, x, y, dimmz, dimmx);

        const real sty = ((C0 * ( sy_current - sy_back1 ) +
                            C1 * ( sy_front1  - sy_back2 ) +
                            C2 * ( sy_front2  - sy_back3 ) +
                            C3 * ( sy_front3  - sy_back4 )) * dyi );

        if (z < nzf && x < nxf)
        {
            vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
        }
    }
}

__global__
void compute_component_vcell_TR_cuda ( real*    vptr,
                                 const real*    szptr,
                                 const real*    sxptr,
                                 const real*    syptr,
                                 const real*    rho,
                                 const real     dt,
                                 const real     dzi,
                                 const real     dxi,
                                 const real     dyi,
                                 const integer  nz0,
                                 const integer  nzf,
                                 const integer  nx0,
                                 const integer  nxf,
                                 const integer  ny0,
                                 const integer  nyf,
                                 const offset_t SZ,
                                 const offset_t SX,
                                 const offset_t SY,
                                 const integer  dimmz,
                                 const integer  dimmx)
{
    for(integer z = blockIdx.x * blockDim.x + threadIdx.x + nz0;
            z < nzf;
            z += gridDim.x * blockDim.x)
    {
        for(integer x = blockIdx.y * blockDim.y + threadIdx.y + nx0;
                x < nxf;
                x += gridDim.y * blockDim.y)
        {
            for(integer y = ny0;
                    y < nyf;
                    y++)
            {
                const real lrho = rho_TR(rho, z, x, y, dimmz, dimmx);

                const real stx  = stencil_X(SX, sxptr, dxi, z, x, y, dimmz, dimmx);
                const real sty  = stencil_Y(SY, syptr, dyi, z, x, y, dimmz, dimmx);
                const real stz  = stencil_Z(SZ, szptr, dzi, z, x, y, dimmz, dimmx);

                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
}

void compute_component_vcell_TR (      real*    vptr,
                                 const real*    szptr,
                                 const real*    sxptr,
                                 const real*    syptr,
                                 const real*    rho,
                                 const real     dt,
                                 const real     dzi,
                                 const real     dxi,
                                 const real     dyi,
                                 const integer  nz0,
                                 const integer  nzf,
                                 const integer  nx0,
                                 const integer  nxf,
                                 const integer  ny0,
                                 const integer  nyf,
                                 const offset_t SZ,
                                 const offset_t SX,
                                 const offset_t SY,
                                 const integer  dimmz,
                                 const integer  dimmx)
{
#if defined(_OPENMP)
    #pragma omp parallel for
#endif /* end pragma _OPENACC */
    for(integer y=ny0; y < nyf; y++)
    {
        for(integer x=nx0; x < nxf; x++)
        {
            for(integer z=nz0; z < nzf; z++)
            {
                const real lrho = rho_TR(rho, z, x, y, dimmz, dimmx);

                const real stx  = stencil_X(SX, sxptr, dxi, z, x, y, dimmz, dimmx);
                const real sty  = stencil_Y(SY, syptr, dyi, z, x, y, dimmz, dimmx);
                const real stz  = stencil_Z(SZ, szptr, dzi, z, x, y, dimmz, dimmx);

                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
};

//////////////////////////////////// VCELL BR /////////////////////////////////////

template <const int HALO  = 4,
          const int BDIMX = 32,
          const int BDIMY = 4,
          const int NX    = BDIMX+2*HALO,
          const int NY    = BDIMY+2*HALO>
__global__
__launch_bounds__(128, 8)
void compute_component_vcell_BR_cuda ( real* __restrict__ vptr,
                                 const real* __restrict__ szptr,
                                 const real* __restrict__ sxptr,
                                 const real* __restrict__ syptr,
                                 const real* __restrict__ rho,
                                 const real               dt,
                                 const real               dzi,
                                 const real               dxi,
                                 const real               dyi,
                                 const integer            nz0,
                                 const integer            nzf,
                                 const integer            nx0,
                                 const integer            nxf,
                                 const integer            ny0,
                                 const integer            nyf,
                                 const offset_t           SZ,
                                 const offset_t           SX,
                                 const offset_t           SY,
                                 const integer            dimmz,
                                 const integer            dimmx)
{
    integer z = blockIdx.x * blockDim.x + threadIdx.x + nz0;
    integer x = blockIdx.y * blockDim.y + threadIdx.y + nx0;

    __shared__ real sx_smem[NY][NX];
    real sy_front1, sy_front2, sy_front3;
    real sy_back1, sy_back2, sy_back3, sy_back4;
    real sy_current;

    sy_back3   = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+0+SY,dimmz,dimmx)] : 0.0f;
    sy_back2   = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+1+SY,dimmz,dimmx)] : 0.0f;
    sy_back1   = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+2+SY,dimmz,dimmx)] : 0.0f;
    sy_current = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+3+SY,dimmz,dimmx)] : 0.0f;
    sy_front1  = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+4+SY,dimmz,dimmx)] : 0.0f;
    sy_front2  = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+5+SY,dimmz,dimmx)] : 0.0f;
    sy_front3  = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+6+SY,dimmz,dimmx)] : 0.0f;

    __shared__ real rho_smem[BDIMY+1][BDIMX+1];
    real rho_current, rho_front1;
    rho_front1 = (z<dimmz && x<dimmx) ? rho[IDX(z,x,ny0,dimmz,dimmx)] : 0.0f;

    for(integer y = ny0;
                y < nyf;
                y++)
    {
        /////////// register tiling-advance plane ////////////////
        sy_back4   = sy_back3;
        sy_back3   = sy_back2;
        sy_back2   = sy_back1;
        sy_back1   = sy_current;
        sy_current = sy_front1;
        sy_front1  = sy_front2;
        sy_front2  = sy_front3;
        sy_front3  = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,y+SY+HALO-1,dimmz,dimmx)] : 0.0f;
        ///////////////////////
        rho_current = rho_front1;
        rho_front1  = (z<dimmz && x<dimmx) ? rho[IDX(z,x,y+1,dimmz,dimmx)] : 0.0f;
        //////////////////////////////////////////////////////////

        const real lrho = rho_BR_smem <BDIMX,BDIMY> (rho_smem, rho, rho_current, rho_front1, z, x, y, dimmz, dimmx);

        const real stz = stencil_Z_shfl <HALO,BDIMX> (SZ, szptr, dzi, z, x, y, dimmz, dimmx);

        const real stx = stencil_X_smem <NX,NY,HALO,BDIMY> (SX, sx_smem, sxptr, dxi, z, x, y, dimmz, dimmx);

        const real sty = ((C0 * ( sy_current - sy_back1 ) +
                            C1 * ( sy_front1  - sy_back2 ) +
                            C2 * ( sy_front2  - sy_back3 ) +
                            C3 * ( sy_front3  - sy_back4 )) * dyi );

        if (z < nzf && x < nxf)
        {
            vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
        }
    }
}


__global__
void compute_component_vcell_BR_cuda ( real*    vptr,
                                 const real*    szptr,
                                 const real*    sxptr,
                                 const real*    syptr,
                                 const real*    rho,
                                 const real     dt,
                                 const real     dzi,
                                 const real     dxi,
                                 const real     dyi,
                                 const integer  nz0,
                                 const integer  nzf,
                                 const integer  nx0,
                                 const integer  nxf,
                                 const integer  ny0,
                                 const integer  nyf,
                                 const offset_t SZ,
                                 const offset_t SX,
                                 const offset_t SY,
                                 const integer  dimmz,
                                 const integer  dimmx)
{
    for(integer z = blockIdx.x * blockDim.x + threadIdx.x + nz0;
            z < nzf;
            z += gridDim.x * blockDim.x)
    {
        for(integer x = blockIdx.y * blockDim.y + threadIdx.y + nx0;
                x < nxf;
                x += gridDim.y * blockDim.y)
        {
            for(integer y = ny0;
                    y < nyf;
                    y++)
            {
                const real lrho = rho_BR(rho, z, x, y, dimmz, dimmx);

                const real stx  = stencil_X(SX, sxptr, dxi, z, x, y, dimmz, dimmx);
                const real sty  = stencil_Y(SY, syptr, dyi, z, x, y, dimmz, dimmx);
                const real stz  = stencil_Z(SZ, szptr, dzi, z, x, y, dimmz, dimmx);

                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
}

void compute_component_vcell_BR (      real*    vptr,
                                 const real*    szptr,
                                 const real*    sxptr,
                                 const real*    syptr,
                                 const real*    rho,
                                 const real     dt,
                                 const real     dzi,
                                 const real     dxi,
                                 const real     dyi,
                                 const integer  nz0,
                                 const integer  nzf,
                                 const integer  nx0,
                                 const integer  nxf,
                                 const integer  ny0,
                                 const integer  nyf,
                                 const offset_t SZ,
                                 const offset_t SX,
                                 const offset_t SY,
                                 const integer  dimmz,
                                 const integer  dimmx)
{
#if defined(_OPENMP)
    #pragma omp parallel for
#endif /* end pragma _OPENACC */
    for(integer y=ny0; y < nyf; y++)
    {
        for(integer x=nx0; x < nxf; x++)
        {
            for(integer z=nz0; z < nzf; z++)
            {
                const real lrho = rho_BR(rho, z, x, y, dimmz, dimmx);

                const real stx  = stencil_X(SX, sxptr, dxi, z, x, y, dimmz, dimmx);
                const real sty  = stencil_Y(SY, syptr, dyi, z, x, y, dimmz, dimmx);
                const real stz  = stencil_Z(SZ, szptr, dzi, z, x, y, dimmz, dimmx);

                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
};
//////////////////////////////////// VCELL BL /////////////////////////////////////

template <const int HALO  = 4,
          const int BDIMX = 32,
          const int BDIMY = 4,
          const int NX    = BDIMX+2*HALO,
          const int NY    = BDIMY+2*HALO>
__global__
__launch_bounds__(128, 8)
void compute_component_vcell_BL_cuda ( real* __restrict__ vptr,
                                 const real* __restrict__ szptr,
                                 const real* __restrict__ sxptr,
                                 const real* __restrict__ syptr,
                                 const real* __restrict__ rho,
                                 const real               dt,
                                 const real               dzi,
                                 const real               dxi,
                                 const real               dyi,
                                 const integer            nz0,
                                 const integer            nzf,
                                 const integer            nx0,
                                 const integer            nxf,
                                 const integer            ny0,
                                 const integer            nyf,
                                 const offset_t           SZ,
                                 const offset_t           SX,
                                 const offset_t           SY,
                                 const integer            dimmz,
                                 const integer            dimmx)
{
    integer z = blockIdx.x * blockDim.x + threadIdx.x + nz0;
    integer x = blockIdx.y * blockDim.y + threadIdx.y + nx0;

    __shared__ real sx_smem[NY][NX];
    real sy_front1, sy_front2, sy_front3;
    real sy_back1, sy_back2, sy_back3, sy_back4;
    real sy_current;

    sy_back3   = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+0+SY,dimmz,dimmx)] : 0.0f;
    sy_back2   = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+1+SY,dimmz,dimmx)] : 0.0f;
    sy_back1   = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+2+SY,dimmz,dimmx)] : 0.0f;
    sy_current = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+3+SY,dimmz,dimmx)] : 0.0f;
    sy_front1  = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+4+SY,dimmz,dimmx)] : 0.0f;
    sy_front2  = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+5+SY,dimmz,dimmx)] : 0.0f;
    sy_front3  = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,ny0-HALO+6+SY,dimmz,dimmx)] : 0.0f;

    for(integer y = ny0;
            y < nyf;
            y++)
    {
        /////////// register tiling-advance plane ////////////////
        sy_back4   = sy_back3;
        sy_back3   = sy_back2;
        sy_back2   = sy_back1;
        sy_back1   = sy_current;
        sy_current = sy_front1;
        sy_front1  = sy_front2;
        sy_front2  = sy_front3;
        sy_front3  = (z<dimmz && x<dimmx) ? syptr[IDX(z,x,y+SY+HALO-1,dimmz,dimmx)] : 0.0f;
        //////////////////////////////////////////////////////////

        const real lrho = rho_BL_shfl <BDIMX> (rho, z, x, y, dimmz, dimmx);

        const real stz = stencil_Z_shfl <HALO,BDIMX> (SZ, szptr, dzi, z, x, y, dimmz, dimmx);

        const real stx = stencil_X_smem <NX,NY,HALO,BDIMY> (SX, sx_smem, sxptr, dxi, z, x, y, dimmz, dimmx);

        const real sty = ((C0 * ( sy_current - sy_back1 ) +
                            C1 * ( sy_front1  - sy_back2 ) +
                            C2 * ( sy_front2  - sy_back3 ) +
                            C3 * ( sy_front3  - sy_back4 )) * dyi );

        if (z < nzf && x < nxf)
        {
            vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
        }
    }
}


__global__
void compute_component_vcell_BL_cuda ( real*    vptr,
                                 const real*    szptr,
                                 const real*    sxptr,
                                 const real*    syptr,
                                 const real*    rho,
                                 const real     dt,
                                 const real     dzi,
                                 const real     dxi,
                                 const real     dyi,
                                 const integer  nz0,
                                 const integer  nzf,
                                 const integer  nx0,
                                 const integer  nxf,
                                 const integer  ny0,
                                 const integer  nyf,
                                 const offset_t SZ,
                                 const offset_t SX,
                                 const offset_t SY,
                                 const integer  dimmz,
                                 const integer  dimmx)
{
    for(integer z = blockIdx.x * blockDim.x + threadIdx.x + nz0;
            z < nzf;
            z += gridDim.x * blockDim.x)
    {
        for(integer x = blockIdx.y * blockDim.y + threadIdx.y + nx0;
                x < nxf;
                x += gridDim.y * blockDim.y)
        {
            for(integer y = ny0;
                    y < nyf;
                    y++)
            {
                const real lrho = rho_BL(rho, z, x, y, dimmz, dimmx);

                const real stx  = stencil_X(SX, sxptr, dxi, z, x, y, dimmz, dimmx);
                const real sty  = stencil_Y(SY, syptr, dyi, z, x, y, dimmz, dimmx);
                const real stz  = stencil_Z(SZ, szptr, dzi, z, x, y, dimmz, dimmx);

                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
}

void compute_component_vcell_BL (      real*    vptr,
                                 const real*    szptr,
                                 const real*    sxptr,
                                 const real*    syptr,
                                 const real*    rho,
                                 const real     dt,
                                 const real     dzi,
                                 const real     dxi,
                                 const real     dyi,
                                 const integer  nz0,
                                 const integer  nzf,
                                 const integer  nx0,
                                 const integer  nxf,
                                 const integer  ny0,
                                 const integer  nyf,
                                 const offset_t SZ,
                                 const offset_t SX,
                                 const offset_t SY,
                                 const integer  dimmz,
                                 const integer  dimmx)
{
#if defined(_OPENMP)
    #pragma omp parallel for
#endif /* end pragma _OPENACC */
    for(integer y=ny0; y < nyf; y++)
    {
        for(integer x=nx0; x < nxf; x++)
        {
            for(integer z=nz0; z < nzf; z++)
            {
                const real lrho = rho_BL(rho, z, x, y, dimmz, dimmx);

                const real stx  = stencil_X(SX, sxptr, dxi, z, x, y, dimmz, dimmx);
                const real sty  = stencil_Y(SY, syptr, dyi, z, x, y, dimmz, dimmx);
                const real stz  = stencil_Z(SZ, szptr, dzi, z, x, y, dimmz, dimmx);

                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
};

void velocity_propagator(v_t           v,
                         s_t           s,
                         coeff_t       coeffs,
                         real*         rho,
                         const real    dt,
                         const real    dzi,
                         const real    dxi,
                         const real    dyi,
                         const integer nz0,
                         const integer nzf,
                         const integer nx0,
                         const integer nxf,
                         const integer ny0,
                         const integer nyf,
                         const integer dimmz,
                         const integer dimmx,
                         const cudaStream_t BR,
                         const cudaStream_t BL,
                         const cudaStream_t TR,
                         const cudaStream_t TL)
{
#if defined(DEBUG)
    fprintf(stderr, "Integration limits of %s are (z "I"-"I",x "I"-"I",y "I"-"I")\n", __FUNCTION__, nz0,nzf,nx0,nxf,ny0,nyf);
#endif

#if !defined(USE_CUDA)

#if defined(__INTEL_COMPILER)
    #pragma forceinline recursive
#endif
    {
        compute_component_vcell_TL (v.tl.w, s.bl.zz, s.tr.xz, s.tl.yz, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, back_offset, forw_offset, dimmz, dimmx);
        compute_component_vcell_TR (v.tr.w, s.br.zz, s.tl.xz, s.tr.yz, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, forw_offset, back_offset, dimmz, dimmx);
        compute_component_vcell_BL (v.bl.w, s.tl.zz, s.br.xz, s.bl.yz, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, back_offset, dimmz, dimmx);
        compute_component_vcell_BR (v.br.w, s.tr.zz, s.bl.xz, s.br.yz, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, forw_offset, forw_offset, dimmz, dimmx);
        compute_component_vcell_TL (v.tl.u, s.bl.xz, s.tr.xx, s.tl.xy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, back_offset, forw_offset, dimmz, dimmx);
        compute_component_vcell_TR (v.tr.u, s.br.xz, s.tl.xx, s.tr.xy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, forw_offset, back_offset, dimmz, dimmx);
        compute_component_vcell_BL (v.bl.u, s.tl.xz, s.br.xx, s.bl.xy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, back_offset, dimmz, dimmx);
        compute_component_vcell_BR (v.br.u, s.tr.xz, s.bl.xx, s.br.xy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, forw_offset, forw_offset, dimmz, dimmx);
        compute_component_vcell_TL (v.tl.v, s.bl.yz, s.tr.xy, s.tl.yy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, back_offset, forw_offset, dimmz, dimmx);
        compute_component_vcell_TR (v.tr.v, s.br.yz, s.tl.xy, s.tr.yy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, forw_offset, back_offset, dimmz, dimmx);
        compute_component_vcell_BL (v.bl.v, s.tl.yz, s.br.xy, s.bl.yy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, back_offset, dimmz, dimmx);
        compute_component_vcell_BR (v.br.v, s.tr.yz, s.bl.xy, s.br.yy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, forw_offset, forw_offset, dimmz, dimmx);
    }

#else /* USE_CUDA -> CALL GPU kernels*/
    const int block_dim_x = 32;
    const int block_dim_y = 4;

    dim3 grid_dim( ((nzf-nz0) + block_dim_x-1)/block_dim_x,
                   ((nxf-nx0) + block_dim_y-1)/block_dim_y,
                      1 );
    dim3 block_dim(block_dim_x, block_dim_y, 1);

    {
#ifdef NAIVE_KERNELS
        compute_component_vcell_TL_cuda <<<grid_dim, block_dim, 0, TL>>>(v.tl.w, s.bl.zz, s.tr.xz, s.tl.yz, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, back_offset, forw_offset, dimmz, dimmx);
        compute_component_vcell_TR_cuda <<<grid_dim, block_dim, 0, TR>>>(v.tr.w, s.br.zz, s.tl.xz, s.tr.yz, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, forw_offset, back_offset, dimmz, dimmx);
        compute_component_vcell_BL_cuda <<<grid_dim, block_dim, 0, BL>>>(v.bl.w, s.tl.zz, s.br.xz, s.bl.yz, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, back_offset, dimmz, dimmx);
        compute_component_vcell_BR_cuda <<<grid_dim, block_dim, 0, BR>>>(v.br.w, s.tr.zz, s.bl.xz, s.br.yz, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, forw_offset, forw_offset, dimmz, dimmx);
        compute_component_vcell_TL_cuda <<<grid_dim, block_dim, 0, TL>>>(v.tl.u, s.bl.xz, s.tr.xx, s.tl.xy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, back_offset, forw_offset, dimmz, dimmx);
        compute_component_vcell_TR_cuda <<<grid_dim, block_dim, 0, TR>>>(v.tr.u, s.br.xz, s.tl.xx, s.tr.xy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, forw_offset, back_offset, dimmz, dimmx);
        compute_component_vcell_BL_cuda <<<grid_dim, block_dim, 0, BL>>>(v.bl.u, s.tl.xz, s.br.xx, s.bl.xy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, back_offset, dimmz, dimmx);
        compute_component_vcell_BR_cuda <<<grid_dim, block_dim, 0, BR>>>(v.br.u, s.tr.xz, s.bl.xx, s.br.xy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, forw_offset, forw_offset, dimmz, dimmx);
        compute_component_vcell_TL_cuda <<<grid_dim, block_dim, 0, TL>>>(v.tl.v, s.bl.yz, s.tr.xy, s.tl.yy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, back_offset, forw_offset, dimmz, dimmx);
        compute_component_vcell_TR_cuda <<<grid_dim, block_dim, 0, TR>>>(v.tr.v, s.br.yz, s.tl.xy, s.tr.yy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, forw_offset, back_offset, dimmz, dimmx);
        compute_component_vcell_BL_cuda <<<grid_dim, block_dim, 0, BL>>>(v.bl.v, s.tl.yz, s.br.xy, s.bl.yy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, back_offset, dimmz, dimmx);
        compute_component_vcell_BR_cuda <<<grid_dim, block_dim, 0, BR>>>(v.br.v, s.tr.yz, s.bl.xy, s.br.yy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, forw_offset, forw_offset, dimmz, dimmx);

#else /*OPTIMIZED KERNELS*/
        compute_component_vcell_TL_cuda <4,block_dim_x,block_dim_y><<<grid_dim, block_dim, 0, TL>>>(v.tl.w, s.bl.zz, s.tr.xz, s.tl.yz, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, back_offset, forw_offset, dimmz, dimmx);
        compute_component_vcell_TR_cuda <4,block_dim_x,block_dim_y><<<grid_dim, block_dim, 0, TR>>>(v.tr.w, s.br.zz, s.tl.xz, s.tr.yz, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, forw_offset, back_offset, dimmz, dimmx);
        compute_component_vcell_BL_cuda <4,block_dim_x,block_dim_y><<<grid_dim, block_dim, 0, BL>>>(v.bl.w, s.tl.zz, s.br.xz, s.bl.yz, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, back_offset, dimmz, dimmx);
        compute_component_vcell_BR_cuda <4,block_dim_x,block_dim_y><<<grid_dim, block_dim, 0, BR>>>(v.br.w, s.tr.zz, s.bl.xz, s.br.yz, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, forw_offset, forw_offset, dimmz, dimmx);
        compute_component_vcell_TL_cuda <4,block_dim_x,block_dim_y><<<grid_dim, block_dim, 0, TL>>>(v.tl.u, s.bl.xz, s.tr.xx, s.tl.xy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, back_offset, forw_offset, dimmz, dimmx);
        compute_component_vcell_TR_cuda <4,block_dim_x,block_dim_y><<<grid_dim, block_dim, 0, TR>>>(v.tr.u, s.br.xz, s.tl.xx, s.tr.xy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, forw_offset, back_offset, dimmz, dimmx);
        compute_component_vcell_BL_cuda <4,block_dim_x,block_dim_y><<<grid_dim, block_dim, 0, BL>>>(v.bl.u, s.tl.xz, s.br.xx, s.bl.xy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, back_offset, dimmz, dimmx);
        compute_component_vcell_BR_cuda <4,block_dim_x,block_dim_y><<<grid_dim, block_dim, 0, BR>>>(v.br.u, s.tr.xz, s.bl.xx, s.br.xy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, forw_offset, forw_offset, dimmz, dimmx);
        compute_component_vcell_TL_cuda <4,block_dim_x,block_dim_y><<<grid_dim, block_dim, 0, TL>>>(v.tl.v, s.bl.yz, s.tr.xy, s.tl.yy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, back_offset, forw_offset, dimmz, dimmx);
        compute_component_vcell_TR_cuda <4,block_dim_x,block_dim_y><<<grid_dim, block_dim, 0, TR>>>(v.tr.v, s.br.yz, s.tl.xy, s.tr.yy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, forw_offset, back_offset, dimmz, dimmx);
        compute_component_vcell_BL_cuda <4,block_dim_x,block_dim_y><<<grid_dim, block_dim, 0, BL>>>(v.bl.v, s.tl.yz, s.br.xy, s.bl.yy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, back_offset, dimmz, dimmx);
        compute_component_vcell_BR_cuda <4,block_dim_x,block_dim_y><<<grid_dim, block_dim, 0, BR>>>(v.br.v, s.tr.yz, s.bl.xy, s.br.yy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, forw_offset, forw_offset, dimmz, dimmx);
#endif
    }

    CUDA_CHECK(cudaGetLastError());
#endif
};





/* ------------------------------------------------------------------------------ */
/*                                                                                */
/*                               CALCULO DE TENSIONES                             */
/*                                                                                */
/* ------------------------------------------------------------------------------ */



HOST_DEVICE_INLINE
void stress_update(real*         sptr,
                   const real    c1,
                   const real    c2,
                   const real    c3,
                   const real    c4,
                   const real    c5,
                   const real    c6,
                   const integer z,
                   const integer x,
                   const integer y,
                   const real    dt,
                   const real    u_x,
                   const real    u_y,
                   const real    u_z,
                   const real    v_x,
                   const real    v_y,
                   const real    v_z,
                   const real    w_x,
                   const real    w_y,
                   const real    w_z,
                   const integer dimmz,
                   const integer dimmx)
{
    real accum  = dt * c1 * u_x;
         accum += dt * c2 * v_y;
         accum += dt * c3 * w_z;
         accum += dt * c4 * (w_y + v_z);
         accum += dt * c5 * (w_x + u_z);
         accum += dt * c6 * (v_x + u_y);
    sptr[IDX(z,x,y,dimmz,dimmx)] += accum;
};

HOST_DEVICE_INLINE
real cell_coeff_BR ( const real*   ptr,
                     const integer z,
                     const integer x,
                     const integer y,
                     const integer dimmz,
                     const integer dimmx)
{
    return ( 1.0f / ( 2.5f  *(ptr[IDX(z  , x  ,y,dimmz,dimmx)] +
                              ptr[IDX(z  , x+1,y,dimmz,dimmx)] +
                              ptr[IDX(z+1, x  ,y,dimmz,dimmx)] +
                              ptr[IDX(z+1, x+1,y,dimmz,dimmx)])) );
};

template <const int NX,
          const int NY,
          const int BDIMX,
          const int BDIMY>
DEVICE inline
real cell_coeff_BR_smem (      real    ptr_smem[NY][NX],
                         const real* __restrict__ ptr_gmem,
                         const integer z,
                         const integer x,
                         const integer y,
                         const integer dimmz,
                         const integer dimmx)
{
    const int tx = threadIdx.x;
    const int ty = threadIdx.y;
    ///////////// intra-block communication///////////////////
    ptr_smem[ty][tx] = ptr_gmem[IDX(z,x,y,dimmz,dimmx)];
    if (tx < 1) ptr_smem[ty][tx+BDIMX] = ptr_gmem[IDX(z+BDIMX,x,y,dimmz,dimmx)];
    if (ty < 1) ptr_smem[ty+BDIMY][tx] = ptr_gmem[IDX(z,x+BDIMY,y,dimmz,dimmx)];
    if (tx == 0 && ty == 0)
                ptr_smem[ty+BDIMY][tx+BDIMX] = ptr_gmem[IDX(z+BDIMX,x+BDIMY,y,dimmz,dimmx)];
    __syncthreads();

    return ( 1.0f / ( 2.5f  *(ptr_smem[ty  ][tx  ] +
                              ptr_smem[ty+1][tx  ] +
                              ptr_smem[ty  ][tx+1] +
                              ptr_smem[ty+1][tx+1])) );
};

HOST_DEVICE_INLINE
real cell_coeff_TL ( const real*   ptr,
                     const integer z,
                     const integer x,
                     const integer y,
                     const integer dimmz,
                     const integer dimmx)
{
    return ( 1.0f / (ptr[IDX(z,x,y,dimmz,dimmx)]));
};

HOST_DEVICE_INLINE
real cell_coeff_BL ( const real*   ptr,
                     const integer z,
                     const integer x,
                     const integer y,
                     const integer dimmz,
                     const integer dimmx)
{
    return ( 1.0f / ( 2.5f *(ptr[IDX(z  ,x,y  ,dimmz,dimmx)] +
                             ptr[IDX(z  ,x,y+1,dimmz,dimmx)] +
                             ptr[IDX(z+1,x,y  ,dimmz,dimmx)] +
                             ptr[IDX(z+1,x,y+1,dimmz,dimmx)])) );
};

template <const int BDIMX>
DEVICE inline
real cell_coeff_BL_shfl (const real* __restrict__ ptr_gmem,
                         const integer z,
                         const integer x,
                         const integer y,
                         const integer dimmz,
                         const integer dimmx)
{
    real current = ptr_gmem[IDX(z,x,y,dimmz,dimmx)];
    real right1  = __shfl_down(current, 1);

    /* For threads without neighbors: */
    if (threadIdx.x >= BDIMX-1 /* 1 */) right1 = ptr_gmem[IDX(z+1,x,y,dimmz,dimmx)];

    real current_front = ptr_gmem[IDX(z,x,y+1,dimmz,dimmx)];
    real right1_front  = __shfl_down(current_front, 1);

    /* For threads without neighbors: */
    if (threadIdx.x >= BDIMX-1 /* 1 */) right1_front = ptr_gmem[IDX(z+1,x,y+1,dimmz,dimmx)];

    return ( 1.0f / ( 2.5f *(current       +
                             current_front +
                             right1        +
                             right1_front)) );
};

HOST_DEVICE_INLINE
real cell_coeff_TR ( const real*   ptr,
                     const integer z,
                     const integer x,
                     const integer y,
                     const integer dimmz,
                     const integer dimmx)
{
    return ( 1.0f / ( 2.5f *(ptr[IDX(z  , x  , y  ,dimmz,dimmx)] +
                             ptr[IDX(z  , x+1, y  ,dimmz,dimmx)] +
                             ptr[IDX(z  , x  , y+1,dimmz,dimmx)] +
                             ptr[IDX(z  , x+1, y+1,dimmz,dimmx)])));
};

template <const int NX,
          const int NY,
          const int BDIMY>
__device__ inline
real cell_coeff_TR_smem (      real    ptr_smem[NY][NX],
                         const real* __restrict__ ptr_gmem,
                         const integer z,
                         const integer x,
                         const integer y,
                         const integer dimmz,
                         const integer dimmx)
{
    const int tx = threadIdx.x;
    const int ty = threadIdx.y;
    ///////////// intra-block communication///////////////////
    ptr_smem[ty][tx] = ptr_gmem[IDX(z,x,y,dimmz,dimmx)];
    if (ty < 1) ptr_smem[ty+BDIMY][tx] = ptr_gmem[IDX(z,x+BDIMY,y,dimmz,dimmx)];
    __syncthreads();

    const real current      = ptr_smem[ty  ][tx];
    const real current_down = ptr_smem[ty+1][tx];

    __syncthreads();
    ptr_smem[ty][tx] = ptr_gmem[IDX(z,x,y+1,dimmz,dimmx)];
    if (ty < 1) ptr_smem[ty+BDIMY][tx] = ptr_gmem[IDX(z,x+BDIMY,y+1,dimmz,dimmx)];
    __syncthreads();

    const real front      = ptr_smem[ty  ][tx];
    const real front_down = ptr_smem[ty+1][tx];

    return ( 1.0f / ( 2.5f *(current      +
                             current_down +
                             front        +
                             front_down)) );
};

HOST_DEVICE_INLINE
real cell_coeff_ARTM_BR( const real*   ptr,
                         const integer z,
                         const integer x,
                         const integer y,
                         const integer dimmz,
                         const integer dimmx)
{
    return ((1.0f / ptr[IDX(z  ,x  ,y,dimmz,dimmx )]  +
             1.0f / ptr[IDX(z  ,x+1,y,dimmz,dimmx )]  +
             1.0f / ptr[IDX(z+1,x  ,y,dimmz,dimmx )]  +
             1.0f / ptr[IDX(z+1,x+1,y,dimmz,dimmx )]) * 0.25f);
};

template <const int NX,
          const int NY,
          const int BDIMX,
          const int BDIMY>
__device__ inline
real cell_coeff_ARTM_BR_smem(      real    ptr_smem[NY][NX],
                             const real* __restrict__ ptr_gmem,
                             const integer z,
                             const integer x,
                             const integer y,
                             const integer dimmz,
                             const integer dimmx)
{
    const int tx = threadIdx.x;
    const int ty = threadIdx.y;
    ///////////// intra-block communication///////////////////
    ptr_smem[ty][tx] = ptr_gmem[IDX(z,x,y,dimmz,dimmx)];
    if (tx < 1) ptr_smem[ty][tx+BDIMX] = ptr_gmem[IDX(z+BDIMX,x,y,dimmz,dimmx)];
    if (ty < 1) ptr_smem[ty+BDIMY][tx] = ptr_gmem[IDX(z,x+BDIMY,y,dimmz,dimmx)];
    if (tx == 0 && ty == 0)
                ptr_smem[ty+BDIMY][tx+BDIMX] = ptr_gmem[IDX(z+BDIMX,x+BDIMY,y,dimmz,dimmx)];
    __syncthreads();

    return ((1.0f / ptr_smem[ty  ][tx  ]  +
             1.0f / ptr_smem[ty+1][tx  ]  +
             1.0f / ptr_smem[ty  ][tx+1]  +
             1.0f / ptr_smem[ty+1][tx+1]) * 0.25f);
};

HOST_DEVICE_INLINE
real cell_coeff_ARTM_TL( const real*   ptr,
                         const integer z,
                         const integer x,
                         const integer y,
                         const integer dimmz,
                         const integer dimmx)
{
    return (1.0f / ptr[IDX(z,x,y,dimmz,dimmx)]);
};

HOST_DEVICE_INLINE
real cell_coeff_ARTM_BL( const real*   ptr,
                         const integer z,
                         const integer x,
                         const integer y,
                         const integer dimmz,
                         const integer dimmx)
{
    return ((1.0f / ptr[IDX(z  ,x,y  ,dimmz,dimmx)]  +
             1.0f / ptr[IDX(z  ,x,y+1,dimmz,dimmx)]  +
             1.0f / ptr[IDX(z+1,x,y  ,dimmz,dimmx)]  +
             1.0f / ptr[IDX(z+1,x,y+1,dimmz,dimmx)]) * 0.25f);
};

template <const int BDIMX>
__device__ inline
real cell_coeff_ARTM_BL_shfl( const float* __restrict__ ptr_gmem,
                              const integer z,
                              const integer x,
                              const integer y,
                              const integer dimmz,
                              const integer dimmx)
{
    real current = ((z < dimmz && x < dimmx) ? ptr_gmem[IDX(z,x,y,dimmz,dimmx)] : 1.0f);
    real right1  = __shfl_down(current, 1);

    /* For threads without neighbors: */
    if (threadIdx.x >= BDIMX-1 /* 1 */) right1 = ((z+1 < dimmz && x < dimmx) ? ptr_gmem[IDX(z+1,x,y,dimmz,dimmx)] : 1.0f);

    real current_front = ((z < dimmz && x < dimmx) ? ptr_gmem[IDX(z,x,y+1,dimmz,dimmx)] : 1.0f);
    real right1_front  = __shfl_down(current_front, 1);

    /* For threads without neighbors: */
    if (threadIdx.x >= BDIMX-1 /* 1 */) right1_front = ((z+1 < dimmz && x < dimmx) ? ptr_gmem[IDX(z+1,x,y+1,dimmz,dimmx)] : 1.0f);


    return ((1.0f / current       +
             1.0f / current_front +
             1.0f / right1        +
             1.0f / right1_front   ) * 0.25f);
};

HOST_DEVICE_INLINE
real cell_coeff_ARTM_TR( const real*   ptr,
                         const integer z,
                         const integer x,
                         const integer y,
                         const integer dimmz,
                         const integer dimmx)
{
    return ((1.0f / ptr[IDX(z,x  ,y  ,dimmz,dimmx)]  +
             1.0f / ptr[IDX(z,x+1,y  ,dimmz,dimmx)]  +
             1.0f / ptr[IDX(z,x  ,y+1,dimmz,dimmx)]  +
             1.0f / ptr[IDX(z,x+1,y+1,dimmz,dimmx)]) * 0.25f);
};

template <const int NX,
          const int NY,
          const int BDIMY>
__device__ inline
real cell_coeff_ARTM_TR_smem(       real    ptr_smem[NY][NX],
                              const real* __restrict__ ptr_gmem,
                              const integer z,
                              const integer x,
                              const integer y,
                              const integer dimmz,
                              const integer dimmx)
{
    const int tx = threadIdx.x;
    const int ty = threadIdx.y;
    ///////////// intra-block communication///////////////////
    ptr_smem[ty][tx] = ptr_gmem[IDX(z,x,y,dimmz,dimmx)];
    if (ty < 1) ptr_smem[ty+BDIMY][tx] = ptr_gmem[IDX(z,x+BDIMY,y,dimmz,dimmx)];
    __syncthreads();

    const real current      = ptr_smem[ty  ][tx];
    const real current_down = ptr_smem[ty+1][tx];

    __syncthreads();
    ptr_smem[ty][tx] = ptr_gmem[IDX(z,x,y+1,dimmz,dimmx)];
    if (ty < 1) ptr_smem[ty+BDIMY][tx] = ptr_gmem[IDX(z,x+BDIMY,y+1,dimmz,dimmx)];
    __syncthreads();

    const real front      = ptr_smem[ty  ][tx];
    const real front_down = ptr_smem[ty+1][tx];

    return ((1.0f / current      +
             1.0f / current_down +
             1.0f / front        +
             1.0f / front_down ) * 0.25f);
};

template <const int HALO  = 4,
          const int BDIMX = 32,
          const int BDIMY =  4,
          const int NX    = BDIMX+2*HALO,
          const int NY    = BDIMY+2*HALO>
__global__
//__launch_bounds__(128, 8)
void compute_component_scell_TR_cuda ( s_t      s,
                                 point_v_t      vnode_z,
                                 point_v_t      vnode_x,
                                 point_v_t      vnode_y,
                                 coeff_t        coeffs,
                                 const real     dt,
                                 const real     dzi,
                                 const real     dxi,
                                 const real     dyi,
                                 const integer  nz0,
                                 const integer  nzf,
                                 const integer  nx0,
                                 const integer  nxf,
                                 const integer  ny0,
                                 const integer  nyf,
                                 const offset_t SZ,
                                 const offset_t SX,
                                 const offset_t SY,
                                 const integer  dimmz,
                                 const integer  dimmx)
{
    const int z = blockIdx.x * blockDim.x + threadIdx.x + nz0;
    const int x = blockIdx.y * blockDim.y + threadIdx.y + nx0;

    __shared__ real bsmem[NY][NX];

    for(int y = ny0; y < nyf; y++)
    {
        real c11, c12, c13, c14, c15, c16;
        real c22, c23, c24, c25, c26;
        real c33, c34, c35, c36;
        real c44, c45, c46;
        real c55, c56;
        real c66;

        if (z < dimmz && x < dimmx)
        {
            c11 = cell_coeff_TR_smem      <NX,NY,BDIMY>(bsmem, coeffs.c11, z, x, y, dimmz, dimmx);
            c12 = cell_coeff_TR_smem      <NX,NY,BDIMY>(bsmem, coeffs.c12, z, x, y, dimmz, dimmx);
            c13 = cell_coeff_TR_smem      <NX,NY,BDIMY>(bsmem, coeffs.c13, z, x, y, dimmz, dimmx);
            c14 = cell_coeff_ARTM_TR_smem <NX,NY,BDIMY>(bsmem, coeffs.c14, z, x, y, dimmz, dimmx);
            c15 = cell_coeff_ARTM_TR_smem <NX,NY,BDIMY>(bsmem, coeffs.c15, z, x, y, dimmz, dimmx);
            c16 = cell_coeff_ARTM_TR_smem <NX,NY,BDIMY>(bsmem, coeffs.c16, z, x, y, dimmz, dimmx);
            c22 = cell_coeff_TR_smem      <NX,NY,BDIMY>(bsmem, coeffs.c22, z, x, y, dimmz, dimmx);
            c23 = cell_coeff_TR_smem      <NX,NY,BDIMY>(bsmem, coeffs.c23, z, x, y, dimmz, dimmx);
            c24 = cell_coeff_ARTM_TR_smem <NX,NY,BDIMY>(bsmem, coeffs.c24, z, x, y, dimmz, dimmx);
            c25 = cell_coeff_ARTM_TR_smem <NX,NY,BDIMY>(bsmem, coeffs.c25, z, x, y, dimmz, dimmx);
            c26 = cell_coeff_ARTM_TR_smem <NX,NY,BDIMY>(bsmem, coeffs.c26, z, x, y, dimmz, dimmx);
            c33 = cell_coeff_TR_smem      <NX,NY,BDIMY>(bsmem, coeffs.c33, z, x, y, dimmz, dimmx);
            c34 = cell_coeff_ARTM_TR_smem <NX,NY,BDIMY>(bsmem, coeffs.c34, z, x, y, dimmz, dimmx);
            c35 = cell_coeff_ARTM_TR_smem <NX,NY,BDIMY>(bsmem, coeffs.c35, z, x, y, dimmz, dimmx);
            c36 = cell_coeff_ARTM_TR_smem <NX,NY,BDIMY>(bsmem, coeffs.c36, z, x, y, dimmz, dimmx);
            c44 = cell_coeff_TR_smem      <NX,NY,BDIMY>(bsmem, coeffs.c44, z, x, y, dimmz, dimmx);
            c45 = cell_coeff_ARTM_TR_smem <NX,NY,BDIMY>(bsmem, coeffs.c45, z, x, y, dimmz, dimmx);
            c46 = cell_coeff_ARTM_TR_smem <NX,NY,BDIMY>(bsmem, coeffs.c46, z, x, y, dimmz, dimmx);
            c55 = cell_coeff_TR_smem      <NX,NY,BDIMY>(bsmem, coeffs.c55, z, x, y, dimmz, dimmx);
            c56 = cell_coeff_ARTM_TR_smem <NX,NY,BDIMY>(bsmem, coeffs.c56, z, x, y, dimmz, dimmx);
            c66 = cell_coeff_TR_smem      <NX,NY,BDIMY>(bsmem, coeffs.c66, z, x, y, dimmz, dimmx);
        }

        const real u_x = stencil_X_smem <NX,NY,HALO,BDIMY>(SX, bsmem, vnode_x.u, dxi, z, x, y, dimmz, dimmx);
        const real v_x = stencil_X_smem <NX,NY,HALO,BDIMY>(SX, bsmem, vnode_x.v, dxi, z, x, y, dimmz, dimmx);
        const real w_x = stencil_X_smem <NX,NY,HALO,BDIMY>(SX, bsmem, vnode_x.w, dxi, z, x, y, dimmz, dimmx);

        real u_y, v_y, w_y;
        if (z < nzf && x < nxf)
        {
            u_y = stencil_Y (SY, vnode_y.u, dyi, z, x, y, dimmz, dimmx);
            v_y = stencil_Y (SY, vnode_y.v, dyi, z, x, y, dimmz, dimmx);
            w_y = stencil_Y (SY, vnode_y.w, dyi, z, x, y, dimmz, dimmx);
        }

        const real u_z = stencil_Z_shfl <HALO,BDIMX>(SZ, vnode_z.u, dzi, z, x, y, dimmz, dimmx);
        const real v_z = stencil_Z_shfl <HALO,BDIMX>(SZ, vnode_z.v, dzi, z, x, y, dimmz, dimmx);
        const real w_z = stencil_Z_shfl <HALO,BDIMX>(SZ, vnode_z.w, dzi, z, x, y, dimmz, dimmx);

        if (z < nzf && x < nxf)
        {
            stress_update (s.tr.xx,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            stress_update (s.tr.yy,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            stress_update (s.tr.zz,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            stress_update (s.tr.yz,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            stress_update (s.tr.xz,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            stress_update (s.tr.xy,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
        }
    }
}


__global__
void compute_component_scell_TR_cuda ( s_t      s,
                                 point_v_t      vnode_z,
                                 point_v_t      vnode_x,
                                 point_v_t      vnode_y,
                                 coeff_t        coeffs,
                                 const real     dt,
                                 const real     dzi,
                                 const real     dxi,
                                 const real     dyi,
                                 const integer  nz0,
                                 const integer  nzf,
                                 const integer  nx0,
                                 const integer  nxf,
                                 const integer  ny0,
                                 const integer  nyf,
                                 const offset_t SZ,
                                 const offset_t SX,
                                 const offset_t SY,
                                 const integer  dimmz,
                                 const integer  dimmx)
{
   for(int z = blockIdx.x * blockDim.x + threadIdx.x + nz0;
            z < nzf;
            z += gridDim.x * blockDim.x)
    {
        for(int x = blockIdx.y * blockDim.y + threadIdx.y + nx0;
                x < nxf;
                x += gridDim.y * blockDim.y)
        {
            for(int y = ny0;
                    y < nyf;
                    y++)
            {
                const real c11 = cell_coeff_TR      (coeffs.c11, z, x, y, dimmz, dimmx);
                const real c12 = cell_coeff_TR      (coeffs.c12, z, x, y, dimmz, dimmx);
                const real c13 = cell_coeff_TR      (coeffs.c13, z, x, y, dimmz, dimmx);
                const real c14 = cell_coeff_ARTM_TR (coeffs.c14, z, x, y, dimmz, dimmx);
                const real c15 = cell_coeff_ARTM_TR (coeffs.c15, z, x, y, dimmz, dimmx);
                const real c16 = cell_coeff_ARTM_TR (coeffs.c16, z, x, y, dimmz, dimmx);
                const real c22 = cell_coeff_TR      (coeffs.c22, z, x, y, dimmz, dimmx);
                const real c23 = cell_coeff_TR      (coeffs.c23, z, x, y, dimmz, dimmx);
                const real c24 = cell_coeff_ARTM_TR (coeffs.c24, z, x, y, dimmz, dimmx);
                const real c25 = cell_coeff_ARTM_TR (coeffs.c25, z, x, y, dimmz, dimmx);
                const real c26 = cell_coeff_ARTM_TR (coeffs.c26, z, x, y, dimmz, dimmx);
                const real c33 = cell_coeff_TR      (coeffs.c33, z, x, y, dimmz, dimmx);
                const real c34 = cell_coeff_ARTM_TR (coeffs.c34, z, x, y, dimmz, dimmx);
                const real c35 = cell_coeff_ARTM_TR (coeffs.c35, z, x, y, dimmz, dimmx);
                const real c36 = cell_coeff_ARTM_TR (coeffs.c36, z, x, y, dimmz, dimmx);
                const real c44 = cell_coeff_TR      (coeffs.c44, z, x, y, dimmz, dimmx);
                const real c45 = cell_coeff_ARTM_TR (coeffs.c45, z, x, y, dimmz, dimmx);
                const real c46 = cell_coeff_ARTM_TR (coeffs.c46, z, x, y, dimmz, dimmx);
                const real c55 = cell_coeff_TR      (coeffs.c55, z, x, y, dimmz, dimmx);
                const real c56 = cell_coeff_ARTM_TR (coeffs.c56, z, x, y, dimmz, dimmx);
                const real c66 = cell_coeff_TR      (coeffs.c66, z, x, y, dimmz, dimmx);

                const real u_x = stencil_X (SX, vnode_x.u, dxi, z, x, y, dimmz, dimmx);
                const real v_x = stencil_X (SX, vnode_x.v, dxi, z, x, y, dimmz, dimmx);
                const real w_x = stencil_X (SX, vnode_x.w, dxi, z, x, y, dimmz, dimmx);

                const real u_y = stencil_Y (SY, vnode_y.u, dyi, z, x, y, dimmz, dimmx);
                const real v_y = stencil_Y (SY, vnode_y.v, dyi, z, x, y, dimmz, dimmx);
                const real w_y = stencil_Y (SY, vnode_y.w, dyi, z, x, y, dimmz, dimmx);

                const real u_z = stencil_Z (SZ, vnode_z.u, dzi, z, x, y, dimmz, dimmx);
                const real v_z = stencil_Z (SZ, vnode_z.v, dzi, z, x, y, dimmz, dimmx);
                const real w_z = stencil_Z (SZ, vnode_z.w, dzi, z, x, y, dimmz, dimmx);

                stress_update (s.tr.xx,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.tr.yy,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.tr.zz,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.tr.yz,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.tr.xz,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.tr.xy,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            }
        }
    }
};

void compute_component_scell_TR (s_t            s,
                                 point_v_t      vnode_z,
                                 point_v_t      vnode_x,
                                 point_v_t      vnode_y,
                                 coeff_t        coeffs,
                                 const real     dt,
                                 const real     dzi,
                                 const real     dxi,
                                 const real     dyi,
                                 const integer  nz0,
                                 const integer  nzf,
                                 const integer  nx0,
                                 const integer  nxf,
                                 const integer  ny0,
                                 const integer  nyf,
                                 const offset_t SZ,
                                 const offset_t SX,
                                 const offset_t SY,
                                 const integer  dimmz,
                                 const integer  dimmx)
{
#if defined(_OPENMP)
    #pragma omp parallel for
#endif /* end pragma _OPENACC */
    for (integer y = ny0; y < nyf; y++)
    {
        for (integer x = nx0; x < nxf; x++)
        {
            for (integer z = nz0; z < nzf; z++ )
            {
                const real c11 = cell_coeff_TR      (coeffs.c11, z, x, y, dimmz, dimmx);
                const real c12 = cell_coeff_TR      (coeffs.c12, z, x, y, dimmz, dimmx);
                const real c13 = cell_coeff_TR      (coeffs.c13, z, x, y, dimmz, dimmx);
                const real c14 = cell_coeff_ARTM_TR (coeffs.c14, z, x, y, dimmz, dimmx);
                const real c15 = cell_coeff_ARTM_TR (coeffs.c15, z, x, y, dimmz, dimmx);
                const real c16 = cell_coeff_ARTM_TR (coeffs.c16, z, x, y, dimmz, dimmx);
                const real c22 = cell_coeff_TR      (coeffs.c22, z, x, y, dimmz, dimmx);
                const real c23 = cell_coeff_TR      (coeffs.c23, z, x, y, dimmz, dimmx);
                const real c24 = cell_coeff_ARTM_TR (coeffs.c24, z, x, y, dimmz, dimmx);
                const real c25 = cell_coeff_ARTM_TR (coeffs.c25, z, x, y, dimmz, dimmx);
                const real c26 = cell_coeff_ARTM_TR (coeffs.c26, z, x, y, dimmz, dimmx);
                const real c33 = cell_coeff_TR      (coeffs.c33, z, x, y, dimmz, dimmx);
                const real c34 = cell_coeff_ARTM_TR (coeffs.c34, z, x, y, dimmz, dimmx);
                const real c35 = cell_coeff_ARTM_TR (coeffs.c35, z, x, y, dimmz, dimmx);
                const real c36 = cell_coeff_ARTM_TR (coeffs.c36, z, x, y, dimmz, dimmx);
                const real c44 = cell_coeff_TR      (coeffs.c44, z, x, y, dimmz, dimmx);
                const real c45 = cell_coeff_ARTM_TR (coeffs.c45, z, x, y, dimmz, dimmx);
                const real c46 = cell_coeff_ARTM_TR (coeffs.c46, z, x, y, dimmz, dimmx);
                const real c55 = cell_coeff_TR      (coeffs.c55, z, x, y, dimmz, dimmx);
                const real c56 = cell_coeff_ARTM_TR (coeffs.c56, z, x, y, dimmz, dimmx);
                const real c66 = cell_coeff_TR      (coeffs.c66, z, x, y, dimmz, dimmx);

                const real u_x = stencil_X (SX, vnode_x.u, dxi, z, x, y, dimmz, dimmx);
                const real v_x = stencil_X (SX, vnode_x.v, dxi, z, x, y, dimmz, dimmx);
                const real w_x = stencil_X (SX, vnode_x.w, dxi, z, x, y, dimmz, dimmx);

                const real u_y = stencil_Y (SY, vnode_y.u, dyi, z, x, y, dimmz, dimmx);
                const real v_y = stencil_Y (SY, vnode_y.v, dyi, z, x, y, dimmz, dimmx);
                const real w_y = stencil_Y (SY, vnode_y.w, dyi, z, x, y, dimmz, dimmx);

                const real u_z = stencil_Z (SZ, vnode_z.u, dzi, z, x, y, dimmz, dimmx);
                const real v_z = stencil_Z (SZ, vnode_z.v, dzi, z, x, y, dimmz, dimmx);
                const real w_z = stencil_Z (SZ, vnode_z.w, dzi, z, x, y, dimmz, dimmx);

                stress_update (s.tr.xx,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.tr.yy,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.tr.zz,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.tr.yz,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.tr.xz,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.tr.xy,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            }
        }
    }
};

template <const int HALO  = 4,
          const int BDIMX = 32,
          const int BDIMY =  4,
          const int NX    = BDIMX+2*HALO,
          const int NY    = BDIMY+2*HALO>
__global__
//__launch_bounds__(128, 8)
void compute_component_scell_TL_cuda ( s_t      s,
                                 point_v_t      vnode_z,
                                 point_v_t      vnode_x,
                                 point_v_t      vnode_y,
                                 coeff_t        coeffs,
                                 const real     dt,
                                 const real     dzi,
                                 const real     dxi,
                                 const real     dyi,
                                 const integer  nz0,
                                 const integer  nzf,
                                 const integer  nx0,
                                 const integer  nxf,
                                 const integer  ny0,
                                 const integer  nyf,
                                 const offset_t SZ,
                                 const offset_t SX,
                                 const offset_t SY,
                                 const integer  dimmz,
                                 const integer  dimmx)
{
    const int z = blockIdx.x * blockDim.x + threadIdx.x + nz0;
    const int x = blockIdx.y * blockDim.y + threadIdx.y + nx0;

    __shared__ real bsmem[NY][NX];

    for(int y = ny0; y < nyf; y++)
    {
        real c11, c12, c13, c14, c15, c16;
        real c22, c23, c24, c25, c26;
        real c33, c34, c35, c36;
        real c44, c45, c46;
        real c55, c56;
        real c66;

        if (z < dimmz && x < dimmx)
        {
            c11 = cell_coeff_TL      (coeffs.c11, z, x, y, dimmz, dimmx);
            c12 = cell_coeff_TL      (coeffs.c12, z, x, y, dimmz, dimmx);
            c13 = cell_coeff_TL      (coeffs.c13, z, x, y, dimmz, dimmx);
            c14 = cell_coeff_ARTM_TL (coeffs.c14, z, x, y, dimmz, dimmx);
            c15 = cell_coeff_ARTM_TL (coeffs.c15, z, x, y, dimmz, dimmx);
            c16 = cell_coeff_ARTM_TL (coeffs.c16, z, x, y, dimmz, dimmx);
            c22 = cell_coeff_TL      (coeffs.c22, z, x, y, dimmz, dimmx);
            c23 = cell_coeff_TL      (coeffs.c23, z, x, y, dimmz, dimmx);
            c24 = cell_coeff_ARTM_TL (coeffs.c24, z, x, y, dimmz, dimmx);
            c25 = cell_coeff_ARTM_TL (coeffs.c25, z, x, y, dimmz, dimmx);
            c26 = cell_coeff_ARTM_TL (coeffs.c26, z, x, y, dimmz, dimmx);
            c33 = cell_coeff_TL      (coeffs.c33, z, x, y, dimmz, dimmx);
            c34 = cell_coeff_ARTM_TL (coeffs.c34, z, x, y, dimmz, dimmx);
            c35 = cell_coeff_ARTM_TL (coeffs.c35, z, x, y, dimmz, dimmx);
            c36 = cell_coeff_ARTM_TL (coeffs.c36, z, x, y, dimmz, dimmx);
            c44 = cell_coeff_TL      (coeffs.c44, z, x, y, dimmz, dimmx);
            c45 = cell_coeff_ARTM_TL (coeffs.c45, z, x, y, dimmz, dimmx);
            c46 = cell_coeff_ARTM_TL (coeffs.c46, z, x, y, dimmz, dimmx);
            c55 = cell_coeff_TL      (coeffs.c55, z, x, y, dimmz, dimmx);
            c56 = cell_coeff_ARTM_TL (coeffs.c56, z, x, y, dimmz, dimmx);
            c66 = cell_coeff_TL      (coeffs.c66, z, x, y, dimmz, dimmx);
        }

        const real u_x = stencil_X_smem <NX,NY,HALO,BDIMY>(SX, bsmem, vnode_x.u, dxi, z, x, y, dimmz, dimmx);
        const real v_x = stencil_X_smem <NX,NY,HALO,BDIMY>(SX, bsmem, vnode_x.v, dxi, z, x, y, dimmz, dimmx);
        const real w_x = stencil_X_smem <NX,NY,HALO,BDIMY>(SX, bsmem, vnode_x.w, dxi, z, x, y, dimmz, dimmx);

        real u_y, v_y, w_y;
        if (z < nzf && x < nxf)
        {
            u_y = stencil_Y (SY, vnode_y.u, dyi, z, x, y, dimmz, dimmx);
            v_y = stencil_Y (SY, vnode_y.v, dyi, z, x, y, dimmz, dimmx);
            w_y = stencil_Y (SY, vnode_y.w, dyi, z, x, y, dimmz, dimmx);
        }

        const real u_z = stencil_Z_shfl <HALO,BDIMX>(SZ, vnode_z.u, dzi, z, x, y, dimmz, dimmx);
        const real v_z = stencil_Z_shfl <HALO,BDIMX>(SZ, vnode_z.v, dzi, z, x, y, dimmz, dimmx);
        const real w_z = stencil_Z_shfl <HALO,BDIMX>(SZ, vnode_z.w, dzi, z, x, y, dimmz, dimmx);

        if (z < nzf && x < nxf)
        {
            stress_update (s.tl.xx,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            stress_update (s.tl.yy,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            stress_update (s.tl.zz,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            stress_update (s.tl.yz,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            stress_update (s.tl.xz,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            stress_update (s.tl.xy,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
        }
    }
}

__global__
void compute_component_scell_TL_cuda ( s_t      s,
                                 point_v_t      vnode_z,
                                 point_v_t      vnode_x,
                                 point_v_t      vnode_y,
                                 coeff_t        coeffs,
                                 const real     dt,
                                 const real     dzi,
                                 const real     dxi,
                                 const real     dyi,
                                 const integer  nz0,
                                 const integer  nzf,
                                 const integer  nx0,
                                 const integer  nxf,
                                 const integer  ny0,
                                 const integer  nyf,
                                 const offset_t SZ,
                                 const offset_t SX,
                                 const offset_t SY,
                                 const integer  dimmz,
                                 const integer  dimmx)
{
   for(int z = blockIdx.x * blockDim.x + threadIdx.x + nz0;
            z < nzf;
            z += gridDim.x * blockDim.x)
    {
        for(int x = blockIdx.y * blockDim.y + threadIdx.y + nx0;
                x < nxf;
                x += gridDim.y * blockDim.y)
        {
            for(int y = ny0;
                    y < nyf;
                    y++)
            {
                const real c11 = cell_coeff_TL      (coeffs.c11, z, x, y, dimmz, dimmx);
                const real c12 = cell_coeff_TL      (coeffs.c12, z, x, y, dimmz, dimmx);
                const real c13 = cell_coeff_TL      (coeffs.c13, z, x, y, dimmz, dimmx);
                const real c14 = cell_coeff_ARTM_TL (coeffs.c14, z, x, y, dimmz, dimmx);
                const real c15 = cell_coeff_ARTM_TL (coeffs.c15, z, x, y, dimmz, dimmx);
                const real c16 = cell_coeff_ARTM_TL (coeffs.c16, z, x, y, dimmz, dimmx);
                const real c22 = cell_coeff_TL      (coeffs.c22, z, x, y, dimmz, dimmx);
                const real c23 = cell_coeff_TL      (coeffs.c23, z, x, y, dimmz, dimmx);
                const real c24 = cell_coeff_ARTM_TL (coeffs.c24, z, x, y, dimmz, dimmx);
                const real c25 = cell_coeff_ARTM_TL (coeffs.c25, z, x, y, dimmz, dimmx);
                const real c26 = cell_coeff_ARTM_TL (coeffs.c26, z, x, y, dimmz, dimmx);
                const real c33 = cell_coeff_TL      (coeffs.c33, z, x, y, dimmz, dimmx);
                const real c34 = cell_coeff_ARTM_TL (coeffs.c34, z, x, y, dimmz, dimmx);
                const real c35 = cell_coeff_ARTM_TL (coeffs.c35, z, x, y, dimmz, dimmx);
                const real c36 = cell_coeff_ARTM_TL (coeffs.c36, z, x, y, dimmz, dimmx);
                const real c44 = cell_coeff_TL      (coeffs.c44, z, x, y, dimmz, dimmx);
                const real c45 = cell_coeff_ARTM_TL (coeffs.c45, z, x, y, dimmz, dimmx);
                const real c46 = cell_coeff_ARTM_TL (coeffs.c46, z, x, y, dimmz, dimmx);
                const real c55 = cell_coeff_TL      (coeffs.c55, z, x, y, dimmz, dimmx);
                const real c56 = cell_coeff_ARTM_TL (coeffs.c56, z, x, y, dimmz, dimmx);
                const real c66 = cell_coeff_TL      (coeffs.c66, z, x, y, dimmz, dimmx);

                const real u_x = stencil_X (SX, vnode_x.u, dxi, z, x, y, dimmz, dimmx);
                const real v_x = stencil_X (SX, vnode_x.v, dxi, z, x, y, dimmz, dimmx);
                const real w_x = stencil_X (SX, vnode_x.w, dxi, z, x, y, dimmz, dimmx);

                const real u_y = stencil_Y (SY, vnode_y.u, dyi, z, x, y, dimmz, dimmx);
                const real v_y = stencil_Y (SY, vnode_y.v, dyi, z, x, y, dimmz, dimmx);
                const real w_y = stencil_Y (SY, vnode_y.w, dyi, z, x, y, dimmz, dimmx);

                const real u_z = stencil_Z (SZ, vnode_z.u, dzi, z, x, y, dimmz, dimmx);
                const real v_z = stencil_Z (SZ, vnode_z.v, dzi, z, x, y, dimmz, dimmx);
                const real w_z = stencil_Z (SZ, vnode_z.w, dzi, z, x, y, dimmz, dimmx);

                stress_update (s.tl.xx,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.tl.yy,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.tl.zz,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.tl.yz,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.tl.xz,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.tl.xy,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            }
        }
    }
};

void compute_component_scell_TL (s_t            s,
                                 point_v_t      vnode_z,
                                 point_v_t      vnode_x,
                                 point_v_t      vnode_y,
                                 coeff_t        coeffs,
                                 const real     dt,
                                 const real     dzi,
                                 const real     dxi,
                                 const real     dyi,
                                 const integer  nz0,
                                 const integer  nzf,
                                 const integer  nx0,
                                 const integer  nxf,
                                 const integer  ny0,
                                 const integer  nyf,
                                 const offset_t SZ,
                                 const offset_t SX,
                                 const offset_t SY,
                                 const integer  dimmz,
                                 const integer  dimmx)
{
#if defined(_OPENMP)
    #pragma omp parallel for
#endif /* end pragma _OPENACC */
    for (integer y = ny0; y < nyf; y++)
    {
        for (integer x = nx0; x < nxf; x++)
        {
            for (integer z = nz0; z < nzf; z++ )
            {
                const real c11 = cell_coeff_TL      (coeffs.c11, z, x, y, dimmz, dimmx);
                const real c12 = cell_coeff_TL      (coeffs.c12, z, x, y, dimmz, dimmx);
                const real c13 = cell_coeff_TL      (coeffs.c13, z, x, y, dimmz, dimmx);
                const real c14 = cell_coeff_ARTM_TL (coeffs.c14, z, x, y, dimmz, dimmx);
                const real c15 = cell_coeff_ARTM_TL (coeffs.c15, z, x, y, dimmz, dimmx);
                const real c16 = cell_coeff_ARTM_TL (coeffs.c16, z, x, y, dimmz, dimmx);
                const real c22 = cell_coeff_TL      (coeffs.c22, z, x, y, dimmz, dimmx);
                const real c23 = cell_coeff_TL      (coeffs.c23, z, x, y, dimmz, dimmx);
                const real c24 = cell_coeff_ARTM_TL (coeffs.c24, z, x, y, dimmz, dimmx);
                const real c25 = cell_coeff_ARTM_TL (coeffs.c25, z, x, y, dimmz, dimmx);
                const real c26 = cell_coeff_ARTM_TL (coeffs.c26, z, x, y, dimmz, dimmx);
                const real c33 = cell_coeff_TL      (coeffs.c33, z, x, y, dimmz, dimmx);
                const real c34 = cell_coeff_ARTM_TL (coeffs.c34, z, x, y, dimmz, dimmx);
                const real c35 = cell_coeff_ARTM_TL (coeffs.c35, z, x, y, dimmz, dimmx);
                const real c36 = cell_coeff_ARTM_TL (coeffs.c36, z, x, y, dimmz, dimmx);
                const real c44 = cell_coeff_TL      (coeffs.c44, z, x, y, dimmz, dimmx);
                const real c45 = cell_coeff_ARTM_TL (coeffs.c45, z, x, y, dimmz, dimmx);
                const real c46 = cell_coeff_ARTM_TL (coeffs.c46, z, x, y, dimmz, dimmx);
                const real c55 = cell_coeff_TL      (coeffs.c55, z, x, y, dimmz, dimmx);
                const real c56 = cell_coeff_ARTM_TL (coeffs.c56, z, x, y, dimmz, dimmx);
                const real c66 = cell_coeff_TL      (coeffs.c66, z, x, y, dimmz, dimmx);

                const real u_x = stencil_X (SX, vnode_x.u, dxi, z, x, y, dimmz, dimmx);
                const real v_x = stencil_X (SX, vnode_x.v, dxi, z, x, y, dimmz, dimmx);
                const real w_x = stencil_X (SX, vnode_x.w, dxi, z, x, y, dimmz, dimmx);

                const real u_y = stencil_Y (SY, vnode_y.u, dyi, z, x, y, dimmz, dimmx);
                const real v_y = stencil_Y (SY, vnode_y.v, dyi, z, x, y, dimmz, dimmx);
                const real w_y = stencil_Y (SY, vnode_y.w, dyi, z, x, y, dimmz, dimmx);

                const real u_z = stencil_Z (SZ, vnode_z.u, dzi, z, x, y, dimmz, dimmx);
                const real v_z = stencil_Z (SZ, vnode_z.v, dzi, z, x, y, dimmz, dimmx);
                const real w_z = stencil_Z (SZ, vnode_z.w, dzi, z, x, y, dimmz, dimmx);

                stress_update (s.tl.xx,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.tl.yy,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.tl.zz,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.tl.yz,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.tl.xz,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.tl.xy,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            }
        }
    }
};



template <const int HALO  = 4,
          const int BDIMX = 32,
          const int BDIMY =  4,
          const int NX    = BDIMX+2*HALO,
          const int NY    = BDIMY+2*HALO>
__global__
//__launch_bounds__(128, 8)
void compute_component_scell_BR_cuda ( s_t      s,
                                 point_v_t      vnode_z,
                                 point_v_t      vnode_x,
                                 point_v_t      vnode_y,
                                 coeff_t        coeffs,
                                 const real     dt,
                                 const real     dzi,
                                 const real     dxi,
                                 const real     dyi,
                                 const integer  nz0,
                                 const integer  nzf,
                                 const integer  nx0,
                                 const integer  nxf,
                                 const integer  ny0,
                                 const integer  nyf,
                                 const offset_t SZ,
                                 const offset_t SX,
                                 const offset_t SY,
                                 const integer  dimmz,
                                 const integer  dimmx)
{
    const int z = blockIdx.x * blockDim.x + threadIdx.x + nz0;
    const int x = blockIdx.y * blockDim.y + threadIdx.y + nx0;

    __shared__ real bsmem[NY][NX];

    for(int y = ny0; y < nyf; y++)
    {
        real c11, c12, c13, c14, c15, c16;
        real c22, c23, c24, c25, c26;
        real c33, c34, c35, c36;
        real c44, c45, c46;
        real c55, c56;
        real c66;

        if (z < dimmz && x < dimmx)
        {
            c11 = cell_coeff_BR_smem      <NX,NY,BDIMX,BDIMY>(bsmem, coeffs.c11, z, x, y, dimmz, dimmx);
            c12 = cell_coeff_BR_smem      <NX,NY,BDIMX,BDIMY>(bsmem, coeffs.c12, z, x, y, dimmz, dimmx);
            c13 = cell_coeff_BR_smem      <NX,NY,BDIMX,BDIMY>(bsmem, coeffs.c13, z, x, y, dimmz, dimmx);
            c14 = cell_coeff_ARTM_BR_smem <NX,NY,BDIMX,BDIMY>(bsmem, coeffs.c14, z, x, y, dimmz, dimmx);
            c15 = cell_coeff_ARTM_BR_smem <NX,NY,BDIMX,BDIMY>(bsmem, coeffs.c15, z, x, y, dimmz, dimmx);
            c16 = cell_coeff_ARTM_BR_smem <NX,NY,BDIMX,BDIMY>(bsmem, coeffs.c16, z, x, y, dimmz, dimmx);
            c22 = cell_coeff_BR_smem      <NX,NY,BDIMX,BDIMY>(bsmem, coeffs.c22, z, x, y, dimmz, dimmx);
            c23 = cell_coeff_BR_smem      <NX,NY,BDIMX,BDIMY>(bsmem, coeffs.c23, z, x, y, dimmz, dimmx);
            c24 = cell_coeff_ARTM_BR_smem <NX,NY,BDIMX,BDIMY>(bsmem, coeffs.c24, z, x, y, dimmz, dimmx);
            c25 = cell_coeff_ARTM_BR_smem <NX,NY,BDIMX,BDIMY>(bsmem, coeffs.c25, z, x, y, dimmz, dimmx);
            c26 = cell_coeff_ARTM_BR_smem <NX,NY,BDIMX,BDIMY>(bsmem, coeffs.c26, z, x, y, dimmz, dimmx);
            c33 = cell_coeff_BR_smem      <NX,NY,BDIMX,BDIMY>(bsmem, coeffs.c33, z, x, y, dimmz, dimmx);
            c34 = cell_coeff_ARTM_BR_smem <NX,NY,BDIMX,BDIMY>(bsmem, coeffs.c34, z, x, y, dimmz, dimmx);
            c35 = cell_coeff_ARTM_BR_smem <NX,NY,BDIMX,BDIMY>(bsmem, coeffs.c35, z, x, y, dimmz, dimmx);
            c36 = cell_coeff_ARTM_BR_smem <NX,NY,BDIMX,BDIMY>(bsmem, coeffs.c36, z, x, y, dimmz, dimmx);
            c44 = cell_coeff_BR_smem      <NX,NY,BDIMX,BDIMY>(bsmem, coeffs.c44, z, x, y, dimmz, dimmx);
            c45 = cell_coeff_ARTM_BR_smem <NX,NY,BDIMX,BDIMY>(bsmem, coeffs.c45, z, x, y, dimmz, dimmx);
            c46 = cell_coeff_ARTM_BR_smem <NX,NY,BDIMX,BDIMY>(bsmem, coeffs.c46, z, x, y, dimmz, dimmx);
            c55 = cell_coeff_BR_smem      <NX,NY,BDIMX,BDIMY>(bsmem, coeffs.c55, z, x, y, dimmz, dimmx);
            c56 = cell_coeff_ARTM_BR_smem <NX,NY,BDIMX,BDIMY>(bsmem, coeffs.c56, z, x, y, dimmz, dimmx);
            c66 = cell_coeff_BR_smem      <NX,NY,BDIMX,BDIMY>(bsmem, coeffs.c66, z, x, y, dimmz, dimmx);
        }

        const real u_x = stencil_X_smem <NX,NY,HALO,BDIMY>(SX, bsmem, vnode_x.u, dxi, z, x, y, dimmz, dimmx);
        const real v_x = stencil_X_smem <NX,NY,HALO,BDIMY>(SX, bsmem, vnode_x.v, dxi, z, x, y, dimmz, dimmx);
        const real w_x = stencil_X_smem <NX,NY,HALO,BDIMY>(SX, bsmem, vnode_x.w, dxi, z, x, y, dimmz, dimmx);

        real u_y, v_y, w_y;
        if (z < nzf && x < nxf)
        {
            u_y = stencil_Y (SY, vnode_y.u, dyi, z, x, y, dimmz, dimmx);
            v_y = stencil_Y (SY, vnode_y.v, dyi, z, x, y, dimmz, dimmx);
            w_y = stencil_Y (SY, vnode_y.w, dyi, z, x, y, dimmz, dimmx);
        }

        const real u_z = stencil_Z_shfl <HALO,BDIMX>(SZ, vnode_z.u, dzi, z, x, y, dimmz, dimmx);
        const real v_z = stencil_Z_shfl <HALO,BDIMX>(SZ, vnode_z.v, dzi, z, x, y, dimmz, dimmx);
        const real w_z = stencil_Z_shfl <HALO,BDIMX>(SZ, vnode_z.w, dzi, z, x, y, dimmz, dimmx);

        if (z < nzf && x < nxf)
        {
            stress_update (s.br.xx,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            stress_update (s.br.yy,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            stress_update (s.br.zz,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            stress_update (s.br.yz,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            stress_update (s.br.xz,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            stress_update (s.br.xy,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
        }
    }
}

__global__
void compute_component_scell_BR_cuda ( s_t      s,
                                 point_v_t      vnode_z,
                                 point_v_t      vnode_x,
                                 point_v_t      vnode_y,
                                 coeff_t        coeffs,
                                 const real     dt,
                                 const real     dzi,
                                 const real     dxi,
                                 const real     dyi,
                                 const integer  nz0,
                                 const integer  nzf,
                                 const integer  nx0,
                                 const integer  nxf,
                                 const integer  ny0,
                                 const integer  nyf,
                                 const offset_t SZ,
                                 const offset_t SX,
                                 const offset_t SY,
                                 const integer  dimmz,
                                 const integer  dimmx)
{
   for(int z = blockIdx.x * blockDim.x + threadIdx.x + nz0;
            z < nzf;
            z += gridDim.x * blockDim.x)
    {
        for(int x = blockIdx.y * blockDim.y + threadIdx.y + nx0;
                x < nxf;
                x += gridDim.y * blockDim.y)
        {
            for(int y = ny0;
                    y < nyf;
                    y++)
            {
                const real c11 = cell_coeff_BR      (coeffs.c11, z, x, y, dimmz, dimmx);
                const real c12 = cell_coeff_BR      (coeffs.c12, z, x, y, dimmz, dimmx);
                const real c13 = cell_coeff_BR      (coeffs.c13, z, x, y, dimmz, dimmx);
                const real c22 = cell_coeff_BR      (coeffs.c22, z, x, y, dimmz, dimmx);
                const real c23 = cell_coeff_BR      (coeffs.c23, z, x, y, dimmz, dimmx);
                const real c33 = cell_coeff_BR      (coeffs.c33, z, x, y, dimmz, dimmx);
                const real c44 = cell_coeff_BR      (coeffs.c44, z, x, y, dimmz, dimmx);
                const real c55 = cell_coeff_BR      (coeffs.c55, z, x, y, dimmz, dimmx);
                const real c66 = cell_coeff_BR      (coeffs.c66, z, x, y, dimmz, dimmx);

                const real c14 = cell_coeff_ARTM_BR (coeffs.c14, z, x, y, dimmz, dimmx);
                const real c15 = cell_coeff_ARTM_BR (coeffs.c15, z, x, y, dimmz, dimmx);
                const real c16 = cell_coeff_ARTM_BR (coeffs.c16, z, x, y, dimmz, dimmx);
                const real c24 = cell_coeff_ARTM_BR (coeffs.c24, z, x, y, dimmz, dimmx);
                const real c25 = cell_coeff_ARTM_BR (coeffs.c25, z, x, y, dimmz, dimmx);
                const real c26 = cell_coeff_ARTM_BR (coeffs.c26, z, x, y, dimmz, dimmx);
                const real c34 = cell_coeff_ARTM_BR (coeffs.c34, z, x, y, dimmz, dimmx);
                const real c35 = cell_coeff_ARTM_BR (coeffs.c35, z, x, y, dimmz, dimmx);
                const real c36 = cell_coeff_ARTM_BR (coeffs.c36, z, x, y, dimmz, dimmx);
                const real c45 = cell_coeff_ARTM_BR (coeffs.c45, z, x, y, dimmz, dimmx);
                const real c46 = cell_coeff_ARTM_BR (coeffs.c46, z, x, y, dimmz, dimmx);
                const real c56 = cell_coeff_ARTM_BR (coeffs.c56, z, x, y, dimmz, dimmx);

                const real u_x = stencil_X (SX, vnode_x.u, dxi, z, x, y, dimmz, dimmx);
                const real v_x = stencil_X (SX, vnode_x.v, dxi, z, x, y, dimmz, dimmx);
                const real w_x = stencil_X (SX, vnode_x.w, dxi, z, x, y, dimmz, dimmx);

                const real u_y = stencil_Y (SY, vnode_y.u, dyi, z, x, y, dimmz, dimmx);
                const real v_y = stencil_Y (SY, vnode_y.v, dyi, z, x, y, dimmz, dimmx);
                const real w_y = stencil_Y (SY, vnode_y.w, dyi, z, x, y, dimmz, dimmx);

                const real u_z = stencil_Z (SZ, vnode_z.u, dzi, z, x, y, dimmz, dimmx);
                const real v_z = stencil_Z (SZ, vnode_z.v, dzi, z, x, y, dimmz, dimmx);
                const real w_z = stencil_Z (SZ, vnode_z.w, dzi, z, x, y, dimmz, dimmx);

                stress_update (s.br.xx,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.br.yy,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.br.zz,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.br.yz,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.br.xz,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.br.xy,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            }
        }
    }
};


void compute_component_scell_BR (s_t            s,
                                 point_v_t      vnode_z,
                                 point_v_t      vnode_x,
                                 point_v_t      vnode_y,
                                 coeff_t        coeffs,
                                 const real     dt,
                                 const real     dzi,
                                 const real     dxi,
                                 const real     dyi,
                                 const integer  nz0,
                                 const integer  nzf,
                                 const integer  nx0,
                                 const integer  nxf,
                                 const integer  ny0,
                                 const integer  nyf,
                                 const offset_t SZ,
                                 const offset_t SX,
                                 const offset_t SY,
                                 const integer  dimmz,
                                 const integer  dimmx)
{
#if defined(_OPENMP)
    #pragma omp parallel for
#endif /* end pragma _OPENACC */
    for (integer y = ny0; y < nyf; y++)
    {
        for (integer x = nx0; x < nxf; x++)
        {
            for (integer z = nz0; z < nzf; z++ )
            {
                const real c11 = cell_coeff_BR      (coeffs.c11, z, x, y, dimmz, dimmx);
                const real c12 = cell_coeff_BR      (coeffs.c12, z, x, y, dimmz, dimmx);
                const real c13 = cell_coeff_BR      (coeffs.c13, z, x, y, dimmz, dimmx);
                const real c22 = cell_coeff_BR      (coeffs.c22, z, x, y, dimmz, dimmx);
                const real c23 = cell_coeff_BR      (coeffs.c23, z, x, y, dimmz, dimmx);
                const real c33 = cell_coeff_BR      (coeffs.c33, z, x, y, dimmz, dimmx);
                const real c44 = cell_coeff_BR      (coeffs.c44, z, x, y, dimmz, dimmx);
                const real c55 = cell_coeff_BR      (coeffs.c55, z, x, y, dimmz, dimmx);
                const real c66 = cell_coeff_BR      (coeffs.c66, z, x, y, dimmz, dimmx);

                const real c14 = cell_coeff_ARTM_BR (coeffs.c14, z, x, y, dimmz, dimmx);
                const real c15 = cell_coeff_ARTM_BR (coeffs.c15, z, x, y, dimmz, dimmx);
                const real c16 = cell_coeff_ARTM_BR (coeffs.c16, z, x, y, dimmz, dimmx);
                const real c24 = cell_coeff_ARTM_BR (coeffs.c24, z, x, y, dimmz, dimmx);
                const real c25 = cell_coeff_ARTM_BR (coeffs.c25, z, x, y, dimmz, dimmx);
                const real c26 = cell_coeff_ARTM_BR (coeffs.c26, z, x, y, dimmz, dimmx);
                const real c34 = cell_coeff_ARTM_BR (coeffs.c34, z, x, y, dimmz, dimmx);
                const real c35 = cell_coeff_ARTM_BR (coeffs.c35, z, x, y, dimmz, dimmx);
                const real c36 = cell_coeff_ARTM_BR (coeffs.c36, z, x, y, dimmz, dimmx);
                const real c45 = cell_coeff_ARTM_BR (coeffs.c45, z, x, y, dimmz, dimmx);
                const real c46 = cell_coeff_ARTM_BR (coeffs.c46, z, x, y, dimmz, dimmx);
                const real c56 = cell_coeff_ARTM_BR (coeffs.c56, z, x, y, dimmz, dimmx);

                const real u_x = stencil_X (SX, vnode_x.u, dxi, z, x, y, dimmz, dimmx);
                const real v_x = stencil_X (SX, vnode_x.v, dxi, z, x, y, dimmz, dimmx);
                const real w_x = stencil_X (SX, vnode_x.w, dxi, z, x, y, dimmz, dimmx);

                const real u_y = stencil_Y (SY, vnode_y.u, dyi, z, x, y, dimmz, dimmx);
                const real v_y = stencil_Y (SY, vnode_y.v, dyi, z, x, y, dimmz, dimmx);
                const real w_y = stencil_Y (SY, vnode_y.w, dyi, z, x, y, dimmz, dimmx);

                const real u_z = stencil_Z (SZ, vnode_z.u, dzi, z, x, y, dimmz, dimmx);
                const real v_z = stencil_Z (SZ, vnode_z.v, dzi, z, x, y, dimmz, dimmx);
                const real w_z = stencil_Z (SZ, vnode_z.w, dzi, z, x, y, dimmz, dimmx);

                stress_update (s.br.xx,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.br.yy,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.br.zz,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.br.yz,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.br.xz,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (s.br.xy,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            }
        }
    }
};

template <const int HALO  = 4,
          const int BDIMX = 32,
          const int BDIMY =  4,
          const int NX    = BDIMX+2*HALO,
          const int NY    = BDIMY+2*HALO>
__global__
//__launch_bounds__(128, 8)
void compute_component_scell_BL_cuda ( s_t      s,
                                 point_v_t      vnode_z,
                                 point_v_t      vnode_x,
                                 point_v_t      vnode_y,
                                 coeff_t        coeffs,
                                 const real     dt,
                                 const real     dzi,
                                 const real     dxi,
                                 const real     dyi,
                                 const integer  nz0,
                                 const integer  nzf,
                                 const integer  nx0,
                                 const integer  nxf,
                                 const integer  ny0,
                                 const integer  nyf,
                                 const offset_t SZ,
                                 const offset_t SX,
                                 const offset_t SY,
                                 const integer  dimmz,
                                 const integer  dimmx)
{
    const int z = blockIdx.x * blockDim.x + threadIdx.x + nz0;
    const int x = blockIdx.y * blockDim.y + threadIdx.y + nx0;

    __shared__ real bsmem[NY][NX];

    for(int y = ny0; y < nyf; y++)
    {
        real c11, c12, c13, c14, c15, c16;
        real c22, c23, c24, c25, c26;
        real c33, c34, c35, c36;
        real c44, c45, c46;
        real c55, c56;
        real c66;

        if (z < dimmz && x < dimmx)
        {
            c11 = cell_coeff_BL_shfl      <BDIMX>(coeffs.c11, z, x, y, dimmz, dimmx);
            c12 = cell_coeff_BL_shfl      <BDIMX>(coeffs.c12, z, x, y, dimmz, dimmx);
            c13 = cell_coeff_BL_shfl      <BDIMX>(coeffs.c13, z, x, y, dimmz, dimmx);
            c14 = cell_coeff_ARTM_BL_shfl <BDIMX>(coeffs.c14, z, x, y, dimmz, dimmx);
            c15 = cell_coeff_ARTM_BL_shfl <BDIMX>(coeffs.c15, z, x, y, dimmz, dimmx);
            c16 = cell_coeff_ARTM_BL_shfl <BDIMX>(coeffs.c16, z, x, y, dimmz, dimmx);
            c22 = cell_coeff_BL_shfl      <BDIMX>(coeffs.c22, z, x, y, dimmz, dimmx);
            c23 = cell_coeff_BL_shfl      <BDIMX>(coeffs.c23, z, x, y, dimmz, dimmx);
            c24 = cell_coeff_ARTM_BL_shfl <BDIMX>(coeffs.c24, z, x, y, dimmz, dimmx);
            c25 = cell_coeff_ARTM_BL_shfl <BDIMX>(coeffs.c25, z, x, y, dimmz, dimmx);
            c26 = cell_coeff_ARTM_BL_shfl <BDIMX>(coeffs.c26, z, x, y, dimmz, dimmx);
            c33 = cell_coeff_BL_shfl      <BDIMX>(coeffs.c33, z, x, y, dimmz, dimmx);
            c34 = cell_coeff_ARTM_BL_shfl <BDIMX>(coeffs.c34, z, x, y, dimmz, dimmx);
            c35 = cell_coeff_ARTM_BL_shfl <BDIMX>(coeffs.c35, z, x, y, dimmz, dimmx);
            c36 = cell_coeff_ARTM_BL_shfl <BDIMX>(coeffs.c36, z, x, y, dimmz, dimmx);
            c44 = cell_coeff_BL_shfl      <BDIMX>(coeffs.c44, z, x, y, dimmz, dimmx);
            c45 = cell_coeff_ARTM_BL_shfl <BDIMX>(coeffs.c45, z, x, y, dimmz, dimmx);
            c46 = cell_coeff_ARTM_BL_shfl <BDIMX>(coeffs.c46, z, x, y, dimmz, dimmx);
            c55 = cell_coeff_BL_shfl      <BDIMX>(coeffs.c55, z, x, y, dimmz, dimmx);
            c56 = cell_coeff_ARTM_BL_shfl <BDIMX>(coeffs.c56, z, x, y, dimmz, dimmx);
            c66 = cell_coeff_BL_shfl      <BDIMX>(coeffs.c66, z, x, y, dimmz, dimmx);
        }

        const real u_x = stencil_X_smem <NX,NY,HALO,BDIMY>(SX, bsmem, vnode_x.u, dxi, z, x, y, dimmz, dimmx);
        const real v_x = stencil_X_smem <NX,NY,HALO,BDIMY>(SX, bsmem, vnode_x.v, dxi, z, x, y, dimmz, dimmx);
        const real w_x = stencil_X_smem <NX,NY,HALO,BDIMY>(SX, bsmem, vnode_x.w, dxi, z, x, y, dimmz, dimmx);

        real u_y, v_y, w_y;
        if (z < nzf && x < nxf)
        {
            u_y = stencil_Y (SY, vnode_y.u, dyi, z, x, y, dimmz, dimmx);
            v_y = stencil_Y (SY, vnode_y.v, dyi, z, x, y, dimmz, dimmx);
            w_y = stencil_Y (SY, vnode_y.w, dyi, z, x, y, dimmz, dimmx);
        }

        const real u_z = stencil_Z_shfl <HALO,BDIMX>(SZ, vnode_z.u, dzi, z, x, y, dimmz, dimmx);
        const real v_z = stencil_Z_shfl <HALO,BDIMX>(SZ, vnode_z.v, dzi, z, x, y, dimmz, dimmx);
        const real w_z = stencil_Z_shfl <HALO,BDIMX>(SZ, vnode_z.w, dzi, z, x, y, dimmz, dimmx);

        if (z < nzf && x < nxf)
        {
            stress_update (s.br.xx,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            stress_update (s.br.yy,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            stress_update (s.br.zz,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            stress_update (s.br.yz,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            stress_update (s.br.xz,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            stress_update (s.br.xy,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
        }
    }
}

__global__
void compute_component_scell_BL_cuda ( s_t      s,
                                 point_v_t      vnode_z,
                                 point_v_t      vnode_x,
                                 point_v_t      vnode_y,
                                 coeff_t        coeffs,
                                 const real     dt,
                                 const real     dzi,
                                 const real     dxi,
                                 const real     dyi,
                                 const integer  nz0,
                                 const integer  nzf,
                                 const integer  nx0,
                                 const integer  nxf,
                                 const integer  ny0,
                                 const integer  nyf,
                                 const offset_t SZ,
                                 const offset_t SX,
                                 const offset_t SY,
                                 const integer  dimmz,
                                 const integer  dimmx)
{
    for(int z = blockIdx.x * blockDim.x + threadIdx.x + nz0;
            z < nzf;
            z += gridDim.x * blockDim.x)
    {
        for(int x = blockIdx.y * blockDim.y + threadIdx.y + nx0;
                x < nxf;
                x += gridDim.y * blockDim.y)
        {
            for(int y = ny0;
                    y < nyf;
                    y++)
            {
                const real c11 = cell_coeff_BL      (coeffs.c11, z, x, y, dimmz, dimmx);
                const real c12 = cell_coeff_BL      (coeffs.c12, z, x, y, dimmz, dimmx);
                const real c13 = cell_coeff_BL      (coeffs.c13, z, x, y, dimmz, dimmx);
                const real c14 = cell_coeff_ARTM_BL (coeffs.c14, z, x, y, dimmz, dimmx);
                const real c15 = cell_coeff_ARTM_BL (coeffs.c15, z, x, y, dimmz, dimmx);
                const real c16 = cell_coeff_ARTM_BL (coeffs.c16, z, x, y, dimmz, dimmx);
                const real c22 = cell_coeff_BL      (coeffs.c22, z, x, y, dimmz, dimmx);
                const real c23 = cell_coeff_BL      (coeffs.c23, z, x, y, dimmz, dimmx);
                const real c24 = cell_coeff_ARTM_BL (coeffs.c24, z, x, y, dimmz, dimmx);
                const real c25 = cell_coeff_ARTM_BL (coeffs.c25, z, x, y, dimmz, dimmx);
                const real c26 = cell_coeff_ARTM_BL (coeffs.c26, z, x, y, dimmz, dimmx);
                const real c33 = cell_coeff_BL      (coeffs.c33, z, x, y, dimmz, dimmx);
                const real c34 = cell_coeff_ARTM_BL (coeffs.c34, z, x, y, dimmz, dimmx);
                const real c35 = cell_coeff_ARTM_BL (coeffs.c35, z, x, y, dimmz, dimmx);
                const real c36 = cell_coeff_ARTM_BL (coeffs.c36, z, x, y, dimmz, dimmx);
                const real c44 = cell_coeff_BL      (coeffs.c44, z, x, y, dimmz, dimmx);
                const real c45 = cell_coeff_ARTM_BL (coeffs.c45, z, x, y, dimmz, dimmx);
                const real c46 = cell_coeff_ARTM_BL (coeffs.c46, z, x, y, dimmz, dimmx);
                const real c55 = cell_coeff_BL      (coeffs.c55, z, x, y, dimmz, dimmx);
                const real c56 = cell_coeff_ARTM_BL (coeffs.c56, z, x, y, dimmz, dimmx);
                const real c66 = cell_coeff_BL      (coeffs.c66, z, x, y, dimmz, dimmx);

                const real u_x = stencil_X (SX, vnode_x.u, dxi, z, x, y, dimmz, dimmx);
                const real v_x = stencil_X (SX, vnode_x.v, dxi, z, x, y, dimmz, dimmx);
                const real w_x = stencil_X (SX, vnode_x.w, dxi, z, x, y, dimmz, dimmx);

                const real u_y = stencil_Y (SY, vnode_y.u, dyi, z, x, y, dimmz, dimmx);
                const real v_y = stencil_Y (SY, vnode_y.v, dyi, z, x, y, dimmz, dimmx);
                const real w_y = stencil_Y (SY, vnode_y.w, dyi, z, x, y, dimmz, dimmx);

                const real u_z = stencil_Z (SZ, vnode_z.u, dzi, z, x, y, dimmz, dimmx);
                const real v_z = stencil_Z (SZ, vnode_z.v, dzi, z, x, y, dimmz, dimmx);
                const real w_z = stencil_Z (SZ, vnode_z.w, dzi, z, x, y, dimmz, dimmx);

                stress_update (s.br.xx,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx);
                stress_update (s.br.yy,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx);
                stress_update (s.br.zz,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx);
                stress_update (s.br.yz,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx);
                stress_update (s.br.xz,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx);
                stress_update (s.br.xy,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx);
            }
        }
    }
};

void compute_component_scell_BL (s_t            s,
                                 point_v_t      vnode_z,
                                 point_v_t      vnode_x,
                                 point_v_t      vnode_y,
                                 coeff_t        coeffs,
                                 const real     dt,
                                 const real     dzi,
                                 const real     dxi,
                                 const real     dyi,
                                 const integer  nz0,
                                 const integer  nzf,
                                 const integer  nx0,
                                 const integer  nxf,
                                 const integer  ny0,
                                 const integer  nyf,
                                 const offset_t SZ,
                                 const offset_t SX,
                                 const offset_t SY,
                                 const integer  dimmz,
                                 const integer  dimmx)
{
#if defined(_OPENMP)
    #pragma omp parallel for
#endif /* end pragma _OPENACC */
    for (integer y = ny0; y < nyf; y++)
    {
        for (integer x = nx0; x < nxf; x++)
        {
            for (integer z = nz0; z < nzf; z++ )
            {
                const real c11 = cell_coeff_BL      (coeffs.c11, z, x, y, dimmz, dimmx);
                const real c12 = cell_coeff_BL      (coeffs.c12, z, x, y, dimmz, dimmx);
                const real c13 = cell_coeff_BL      (coeffs.c13, z, x, y, dimmz, dimmx);
                const real c14 = cell_coeff_ARTM_BL (coeffs.c14, z, x, y, dimmz, dimmx);
                const real c15 = cell_coeff_ARTM_BL (coeffs.c15, z, x, y, dimmz, dimmx);
                const real c16 = cell_coeff_ARTM_BL (coeffs.c16, z, x, y, dimmz, dimmx);
                const real c22 = cell_coeff_BL      (coeffs.c22, z, x, y, dimmz, dimmx);
                const real c23 = cell_coeff_BL      (coeffs.c23, z, x, y, dimmz, dimmx);
                const real c24 = cell_coeff_ARTM_BL (coeffs.c24, z, x, y, dimmz, dimmx);
                const real c25 = cell_coeff_ARTM_BL (coeffs.c25, z, x, y, dimmz, dimmx);
                const real c26 = cell_coeff_ARTM_BL (coeffs.c26, z, x, y, dimmz, dimmx);
                const real c33 = cell_coeff_BL      (coeffs.c33, z, x, y, dimmz, dimmx);
                const real c34 = cell_coeff_ARTM_BL (coeffs.c34, z, x, y, dimmz, dimmx);
                const real c35 = cell_coeff_ARTM_BL (coeffs.c35, z, x, y, dimmz, dimmx);
                const real c36 = cell_coeff_ARTM_BL (coeffs.c36, z, x, y, dimmz, dimmx);
                const real c44 = cell_coeff_BL      (coeffs.c44, z, x, y, dimmz, dimmx);
                const real c45 = cell_coeff_ARTM_BL (coeffs.c45, z, x, y, dimmz, dimmx);
                const real c46 = cell_coeff_ARTM_BL (coeffs.c46, z, x, y, dimmz, dimmx);
                const real c55 = cell_coeff_BL      (coeffs.c55, z, x, y, dimmz, dimmx);
                const real c56 = cell_coeff_ARTM_BL (coeffs.c56, z, x, y, dimmz, dimmx);
                const real c66 = cell_coeff_BL      (coeffs.c66, z, x, y, dimmz, dimmx);

                const real u_x = stencil_X (SX, vnode_x.u, dxi, z, x, y, dimmz, dimmx);
                const real v_x = stencil_X (SX, vnode_x.v, dxi, z, x, y, dimmz, dimmx);
                const real w_x = stencil_X (SX, vnode_x.w, dxi, z, x, y, dimmz, dimmx);

                const real u_y = stencil_Y (SY, vnode_y.u, dyi, z, x, y, dimmz, dimmx);
                const real v_y = stencil_Y (SY, vnode_y.v, dyi, z, x, y, dimmz, dimmx);
                const real w_y = stencil_Y (SY, vnode_y.w, dyi, z, x, y, dimmz, dimmx);

                const real u_z = stencil_Z (SZ, vnode_z.u, dzi, z, x, y, dimmz, dimmx);
                const real v_z = stencil_Z (SZ, vnode_z.v, dzi, z, x, y, dimmz, dimmx);
                const real w_z = stencil_Z (SZ, vnode_z.w, dzi, z, x, y, dimmz, dimmx);

                stress_update (s.br.xx,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx);
                stress_update (s.br.yy,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx);
                stress_update (s.br.zz,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx);
                stress_update (s.br.yz,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx);
                stress_update (s.br.xz,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx);
                stress_update (s.br.xy,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx);
            }
        }
    }
};


void stress_propagator(s_t           s,
                       v_t           v,
                       coeff_t       coeffs,
                       real*         rho,
                       const real    dt,
                       const real    dzi,
                       const real    dxi,
                       const real    dyi,
                       const integer nz0,
                       const integer nzf,
                       const integer nx0,
                       const integer nxf,
                       const integer ny0,
                       const integer nyf,
                       const integer dimmz,
                       const integer dimmx,
                       const cudaStream_t BR,
                       const cudaStream_t BL,
                       const cudaStream_t TR,
                       const cudaStream_t TL)
{
#if !defined(USE_CUDA)

#if defined(__INTEL_COMPILER)
    #pragma forceinline recursive
#endif
    {
        compute_component_scell_BR ( s, v.tr, v.bl, v.br, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, back_offset, dimmz, dimmx);
        compute_component_scell_BL ( s, v.tl, v.br, v.bl, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, forw_offset, dimmz, dimmx);
        compute_component_scell_TR ( s, v.br, v.tl, v.tr, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, forw_offset, forw_offset, dimmz, dimmx);
        compute_component_scell_TL ( s, v.bl, v.tr, v.tl, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, back_offset, back_offset, dimmz, dimmx);
    }

#else /* USE_CUDA -> CALL CUDA KERNELS */
    const int block_dim_x = 32;
    const int block_dim_y = 4;

    dim3 grid_dim( ((nzf-nz0) + block_dim_x-1)/block_dim_x,
                   ((nxf-nx0) + block_dim_y-1)/block_dim_y,
                      1 );
    dim3 block_dim(block_dim_x, block_dim_y, 1);

    {
#if 1
        compute_component_scell_BR_cuda <<<grid_dim,block_dim, 0, BR>>>( s, v.tr, v.bl, v.br, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, back_offset, dimmz, dimmx);
        compute_component_scell_BL_cuda <<<grid_dim,block_dim, 0, BL>>>( s, v.tl, v.br, v.bl, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, forw_offset, dimmz, dimmx);
        compute_component_scell_TR_cuda <<<grid_dim,block_dim, 0, TR>>>( s, v.br, v.tl, v.tr, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, forw_offset, forw_offset, dimmz, dimmx);
        compute_component_scell_TL_cuda <<<grid_dim,block_dim, 0, TL>>>( s, v.bl, v.tr, v.tl, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, back_offset, back_offset, dimmz, dimmx);
#else /* optimized kernels*/
        compute_component_scell_BR_cuda <4,block_dim_x,block_dim_y><<<grid_dim,block_dim, 0, BR>>>( s, v.tr, v.bl, v.br, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, back_offset, dimmz, dimmx);
        compute_component_scell_BL_cuda <4,block_dim_x,block_dim_y><<<grid_dim,block_dim, 0, BL>>>( s, v.tl, v.br, v.bl, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, forw_offset, dimmz, dimmx);
        compute_component_scell_TR_cuda <4,block_dim_x,block_dim_y><<<grid_dim,block_dim, 0, TR>>>( s, v.br, v.tl, v.tr, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, forw_offset, forw_offset, dimmz, dimmx);
        compute_component_scell_TL_cuda <4,block_dim_x,block_dim_y><<<grid_dim,block_dim, 0, TL>>>( s, v.bl, v.tr, v.tl, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, back_offset, back_offset, dimmz, dimmx);

#endif
    }

    CUDA_CHECK(cudaGetLastError());

#endif /* USE_CUDA */
};

