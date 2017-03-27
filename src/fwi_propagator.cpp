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

#if defined(USE_CUDA)
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
#endif

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
                                 const offset_t _SZ,
                                 const offset_t _SX,
                                 const offset_t _SY,
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

                const real stx  = stencil_X( _SX, sxptr, dxi, z, x, y, dimmz, dimmx);
                const real sty  = stencil_Y( _SY, syptr, dyi, z, x, y, dimmz, dimmx);
                const real stz  = stencil_Z( _SZ, szptr, dzi, z, x, y, dimmz, dimmx);

                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
};

#if defined(USE_CUDA)
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

                const real stx  = stencil_X( SX, sxptr, dxi, z, x, y, dimmz, dimmx);
                const real sty  = stencil_Y( SY, syptr, dyi, z, x, y, dimmz, dimmx);
                const real stz  = stencil_Z( SZ, szptr, dzi, z, x, y, dimmz, dimmx);

                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
}
#endif /* USE_CUDA */

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
                                 const offset_t _SZ,
                                 const offset_t _SX,
                                 const offset_t _SY,
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

                const real stx  = stencil_X( _SX, sxptr, dxi, z, x, y, dimmz, dimmx);
                const real sty  = stencil_Y( _SY, syptr, dyi, z, x, y, dimmz, dimmx);
                const real stz  = stencil_Z( _SZ, szptr, dzi, z, x, y, dimmz, dimmx);

                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
};

#if defined(USE_CUDA)
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

                const real stx  = stencil_X( SX, sxptr, dxi, z, x, y, dimmz, dimmx);
                const real sty  = stencil_Y( SY, syptr, dyi, z, x, y, dimmz, dimmx);
                const real stz  = stencil_Z( SZ, szptr, dzi, z, x, y, dimmz, dimmx);

                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
}
#endif /*USE_CUDA*/

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
                                 const offset_t _SZ,
                                 const offset_t _SX,
                                 const offset_t _SY,
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

                const real stx  = stencil_X( _SX, sxptr, dxi, z, x, y, dimmz, dimmx );
                const real sty  = stencil_Y( _SY, syptr, dyi, z, x, y, dimmz, dimmx );
                const real stz  = stencil_Z( _SZ, szptr, dzi, z, x, y, dimmz, dimmx );

                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
};

#if defined(USE_CUDA)
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

                const real stx  = stencil_X( SX, sxptr, dxi, z, x, y, dimmz, dimmx);
                const real sty  = stencil_Y( SY, syptr, dyi, z, x, y, dimmz, dimmx);
                const real stz  = stencil_Z( SZ, szptr, dzi, z, x, y, dimmz, dimmx);

                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
}
#endif /* USE_CUDA */

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
                                 const offset_t _SZ,
                                 const offset_t _SX,
                                 const offset_t _SY,
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

                const real stx  = stencil_X( _SX, sxptr, dxi, z, x, y, dimmz, dimmx);
                const real sty  = stencil_Y( _SY, syptr, dyi, z, x, y, dimmz, dimmx);
                const real stz  = stencil_Z( _SZ, szptr, dzi, z, x, y, dimmz, dimmx);

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

#if defined(USE_CUDA)
__global__
void compute_component_scell_TR_cuda (s_t        s,
                                 point_v_t       vnode_z,
                                 point_v_t       vnode_x,
                                 point_v_t       vnode_y,
                                 coeff_t         coeffs,
                                 const real      dt,
                                 const real      dzi,
                                 const real      dxi,
                                 const real      dyi,
                                 const integer   nz0,
                                 const integer   nzf,
                                 const integer   nx0,
                                 const integer   nxf,
                                 const integer   ny0,
                                 const integer   nyf,
                                 const offset_t _SZ,
                                 const offset_t _SX,
                                 const offset_t _SY,
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

                const real u_x = stencil_X (_SX, vnode_x.u, dxi, z, x, y, dimmz, dimmx);
                const real v_x = stencil_X (_SX, vnode_x.v, dxi, z, x, y, dimmz, dimmx);
                const real w_x = stencil_X (_SX, vnode_x.w, dxi, z, x, y, dimmz, dimmx);

                const real u_y = stencil_Y (_SY, vnode_y.u, dyi, z, x, y, dimmz, dimmx);
                const real v_y = stencil_Y (_SY, vnode_y.v, dyi, z, x, y, dimmz, dimmx);
                const real w_y = stencil_Y (_SY, vnode_y.w, dyi, z, x, y, dimmz, dimmx);

                const real u_z = stencil_Z (_SZ, vnode_z.u, dzi, z, x, y, dimmz, dimmx);
                const real v_z = stencil_Z (_SZ, vnode_z.v, dzi, z, x, y, dimmz, dimmx);
                const real w_z = stencil_Z (_SZ, vnode_z.w, dzi, z, x, y, dimmz, dimmx);

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
#endif /* USE_CUDA */

void compute_component_scell_TR (s_t             s,
                                 point_v_t       vnode_z,
                                 point_v_t       vnode_x,
                                 point_v_t       vnode_y,
                                 coeff_t         coeffs,
                                 const real      dt,
                                 const real      dzi,
                                 const real      dxi,
                                 const real      dyi,
                                 const integer   nz0,
                                 const integer   nzf,
                                 const integer   nx0,
                                 const integer   nxf,
                                 const integer   ny0,
                                 const integer   nyf,
                                 const offset_t _SZ,
                                 const offset_t _SX,
                                 const offset_t _SY,
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

                const real u_x = stencil_X (_SX, vnode_x.u, dxi, z, x, y, dimmz, dimmx);
                const real v_x = stencil_X (_SX, vnode_x.v, dxi, z, x, y, dimmz, dimmx);
                const real w_x = stencil_X (_SX, vnode_x.w, dxi, z, x, y, dimmz, dimmx);

                const real u_y = stencil_Y (_SY, vnode_y.u, dyi, z, x, y, dimmz, dimmx);
                const real v_y = stencil_Y (_SY, vnode_y.v, dyi, z, x, y, dimmz, dimmx);
                const real w_y = stencil_Y (_SY, vnode_y.w, dyi, z, x, y, dimmz, dimmx);

                const real u_z = stencil_Z (_SZ, vnode_z.u, dzi, z, x, y, dimmz, dimmx);
                const real v_z = stencil_Z (_SZ, vnode_z.v, dzi, z, x, y, dimmz, dimmx);
                const real w_z = stencil_Z (_SZ, vnode_z.w, dzi, z, x, y, dimmz, dimmx);

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

#if defined(USE_CUDA)
__global__
void compute_component_scell_TL_cuda ( s_t       s,
                                 point_v_t       vnode_z,
                                 point_v_t       vnode_x,
                                 point_v_t       vnode_y,
                                 coeff_t         coeffs,
                                 const real      dt,
                                 const real      dzi,
                                 const real      dxi,
                                 const real      dyi,
                                 const integer   nz0,
                                 const integer   nzf,
                                 const integer   nx0,
                                 const integer   nxf,
                                 const integer   ny0,
                                 const integer   nyf,
                                 const offset_t _SZ,
                                 const offset_t _SX,
                                 const offset_t _SY,
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

                const real u_x = stencil_X (_SX, vnode_x.u, dxi, z, x, y, dimmz, dimmx);
                const real v_x = stencil_X (_SX, vnode_x.v, dxi, z, x, y, dimmz, dimmx);
                const real w_x = stencil_X (_SX, vnode_x.w, dxi, z, x, y, dimmz, dimmx);

                const real u_y = stencil_Y (_SY, vnode_y.u, dyi, z, x, y, dimmz, dimmx);
                const real v_y = stencil_Y (_SY, vnode_y.v, dyi, z, x, y, dimmz, dimmx);
                const real w_y = stencil_Y (_SY, vnode_y.w, dyi, z, x, y, dimmz, dimmx);

                const real u_z = stencil_Z (_SZ, vnode_z.u, dzi, z, x, y, dimmz, dimmx);
                const real v_z = stencil_Z (_SZ, vnode_z.v, dzi, z, x, y, dimmz, dimmx);
                const real w_z = stencil_Z (_SZ, vnode_z.w, dzi, z, x, y, dimmz, dimmx);

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
#endif /* USE_CUDA */

void compute_component_scell_TL (s_t             s,
                                 point_v_t       vnode_z,
                                 point_v_t       vnode_x,
                                 point_v_t       vnode_y,
                                 coeff_t         coeffs,
                                 const real      dt,
                                 const real      dzi,
                                 const real      dxi,
                                 const real      dyi,
                                 const integer   nz0,
                                 const integer   nzf,
                                 const integer   nx0,
                                 const integer   nxf,
                                 const integer   ny0,
                                 const integer   nyf,
                                 const offset_t _SZ,
                                 const offset_t _SX,
                                 const offset_t _SY,
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

                const real u_x = stencil_X (_SX, vnode_x.u, dxi, z, x, y, dimmz, dimmx);
                const real v_x = stencil_X (_SX, vnode_x.v, dxi, z, x, y, dimmz, dimmx);
                const real w_x = stencil_X (_SX, vnode_x.w, dxi, z, x, y, dimmz, dimmx);

                const real u_y = stencil_Y (_SY, vnode_y.u, dyi, z, x, y, dimmz, dimmx);
                const real v_y = stencil_Y (_SY, vnode_y.v, dyi, z, x, y, dimmz, dimmx);
                const real w_y = stencil_Y (_SY, vnode_y.w, dyi, z, x, y, dimmz, dimmx);

                const real u_z = stencil_Z (_SZ, vnode_z.u, dzi, z, x, y, dimmz, dimmx);
                const real v_z = stencil_Z (_SZ, vnode_z.v, dzi, z, x, y, dimmz, dimmx);
                const real w_z = stencil_Z (_SZ, vnode_z.w, dzi, z, x, y, dimmz, dimmx);

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

#if defined(USE_CUDA)
__global__
void compute_component_scell_BR_cuda ( s_t       s,
                                 point_v_t       vnode_z,
                                 point_v_t       vnode_x,
                                 point_v_t       vnode_y,
                                 coeff_t         coeffs,
                                 const real      dt,
                                 const real      dzi,
                                 const real      dxi,
                                 const real      dyi,
                                 const integer   nz0,
                                 const integer   nzf,
                                 const integer   nx0,
                                 const integer   nxf,
                                 const integer   ny0,
                                 const integer   nyf,
                                 const offset_t _SZ,
                                 const offset_t _SX,
                                 const offset_t _SY,
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

                const real u_x = stencil_X (_SX, vnode_x.u, dxi, z, x, y, dimmz, dimmx);
                const real v_x = stencil_X (_SX, vnode_x.v, dxi, z, x, y, dimmz, dimmx);
                const real w_x = stencil_X (_SX, vnode_x.w, dxi, z, x, y, dimmz, dimmx);

                const real u_y = stencil_Y (_SY, vnode_y.u, dyi, z, x, y, dimmz, dimmx);
                const real v_y = stencil_Y (_SY, vnode_y.v, dyi, z, x, y, dimmz, dimmx);
                const real w_y = stencil_Y (_SY, vnode_y.w, dyi, z, x, y, dimmz, dimmx);

                const real u_z = stencil_Z (_SZ, vnode_z.u, dzi, z, x, y, dimmz, dimmx);
                const real v_z = stencil_Z (_SZ, vnode_z.v, dzi, z, x, y, dimmz, dimmx);
                const real w_z = stencil_Z (_SZ, vnode_z.w, dzi, z, x, y, dimmz, dimmx);

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
#endif /* USE_CUDA */


void compute_component_scell_BR (s_t             s,
                                 point_v_t       vnode_z,
                                 point_v_t       vnode_x,
                                 point_v_t       vnode_y,
                                 coeff_t         coeffs,
                                 const real      dt,
                                 const real      dzi,
                                 const real      dxi,
                                 const real      dyi,
                                 const integer   nz0,
                                 const integer   nzf,
                                 const integer   nx0,
                                 const integer   nxf,
                                 const integer   ny0,
                                 const integer   nyf,
                                 const offset_t _SZ,
                                 const offset_t _SX,
                                 const offset_t _SY,
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

                const real u_x = stencil_X (_SX, vnode_x.u, dxi, z, x, y, dimmz, dimmx);
                const real v_x = stencil_X (_SX, vnode_x.v, dxi, z, x, y, dimmz, dimmx);
                const real w_x = stencil_X (_SX, vnode_x.w, dxi, z, x, y, dimmz, dimmx);

                const real u_y = stencil_Y (_SY, vnode_y.u, dyi, z, x, y, dimmz, dimmx);
                const real v_y = stencil_Y (_SY, vnode_y.v, dyi, z, x, y, dimmz, dimmx);
                const real w_y = stencil_Y (_SY, vnode_y.w, dyi, z, x, y, dimmz, dimmx);

                const real u_z = stencil_Z (_SZ, vnode_z.u, dzi, z, x, y, dimmz, dimmx);
                const real v_z = stencil_Z (_SZ, vnode_z.v, dzi, z, x, y, dimmz, dimmx);
                const real w_z = stencil_Z (_SZ, vnode_z.w, dzi, z, x, y, dimmz, dimmx);

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

#if defined(USE_CUDA)
__global__
void compute_component_scell_BL_cuda ( s_t       s,
                                 point_v_t       vnode_z,
                                 point_v_t       vnode_x,
                                 point_v_t       vnode_y,
                                 coeff_t         coeffs,
                                 const real      dt,
                                 const real      dzi,
                                 const real      dxi,
                                 const real      dyi,
                                 const integer   nz0,
                                 const integer   nzf,
                                 const integer   nx0,
                                 const integer   nxf,
                                 const integer   ny0,
                                 const integer   nyf,
                                 const offset_t _SZ,
                                 const offset_t _SX,
                                 const offset_t _SY,
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

                const real u_x = stencil_X (_SX, vnode_x.u, dxi, z, x, y, dimmz, dimmx);
                const real v_x = stencil_X (_SX, vnode_x.v, dxi, z, x, y, dimmz, dimmx);
                const real w_x = stencil_X (_SX, vnode_x.w, dxi, z, x, y, dimmz, dimmx);

                const real u_y = stencil_Y (_SY, vnode_y.u, dyi, z, x, y, dimmz, dimmx);
                const real v_y = stencil_Y (_SY, vnode_y.v, dyi, z, x, y, dimmz, dimmx);
                const real w_y = stencil_Y (_SY, vnode_y.w, dyi, z, x, y, dimmz, dimmx);

                const real u_z = stencil_Z (_SZ, vnode_z.u, dzi, z, x, y, dimmz, dimmx);
                const real v_z = stencil_Z (_SZ, vnode_z.v, dzi, z, x, y, dimmz, dimmx);
                const real w_z = stencil_Z (_SZ, vnode_z.w, dzi, z, x, y, dimmz, dimmx);

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
#endif /* USE_CUDA */

void compute_component_scell_BL (s_t             s,
                                 point_v_t       vnode_z,
                                 point_v_t       vnode_x,
                                 point_v_t       vnode_y,
                                 coeff_t         coeffs,
                                 const real      dt,
                                 const real      dzi,
                                 const real      dxi,
                                 const real      dyi,
                                 const integer   nz0,
                                 const integer   nzf,
                                 const integer   nx0,
                                 const integer   nxf,
                                 const integer   ny0,
                                 const integer   nyf,
                                 const offset_t _SZ,
                                 const offset_t _SX,
                                 const offset_t _SY,
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

                const real u_x = stencil_X (_SX, vnode_x.u, dxi, z, x, y, dimmz, dimmx);
                const real v_x = stencil_X (_SX, vnode_x.v, dxi, z, x, y, dimmz, dimmx);
                const real w_x = stencil_X (_SX, vnode_x.w, dxi, z, x, y, dimmz, dimmx);

                const real u_y = stencil_Y (_SY, vnode_y.u, dyi, z, x, y, dimmz, dimmx);
                const real v_y = stencil_Y (_SY, vnode_y.v, dyi, z, x, y, dimmz, dimmx);
                const real w_y = stencil_Y (_SY, vnode_y.w, dyi, z, x, y, dimmz, dimmx);

                const real u_z = stencil_Z (_SZ, vnode_z.u, dzi, z, x, y, dimmz, dimmx);
                const real v_z = stencil_Z (_SZ, vnode_z.v, dzi, z, x, y, dimmz, dimmx);
                const real w_z = stencil_Z (_SZ, vnode_z.w, dzi, z, x, y, dimmz, dimmx);

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
        compute_component_scell_BR_cuda <<<grid_dim,block_dim, 0, BR>>>( s, v.tr, v.bl, v.br, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, back_offset, dimmz, dimmx);
        compute_component_scell_BL_cuda <<<grid_dim,block_dim, 0, BL>>>( s, v.tl, v.br, v.bl, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, forw_offset, dimmz, dimmx);
        compute_component_scell_TR_cuda <<<grid_dim,block_dim, 0, TR>>>( s, v.br, v.tl, v.tr, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, forw_offset, forw_offset, dimmz, dimmx);
        compute_component_scell_TL_cuda <<<grid_dim,block_dim, 0, TL>>>( s, v.bl, v.tr, v.tl, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, back_offset, back_offset, dimmz, dimmx);
    }

    CUDA_CHECK(cudaGetLastError());

#endif /* USE_CUDA */
};

