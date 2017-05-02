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

#include "fwi/fwi_propagator.h"

#define IDX(z,x,y,dimmz,dimmx) (((y*dimmx)+x)*dimmz+z)

//inline
//integer IDX (const integer z,
//             const integer x,
//             const integer y,
//             const integer dimmz,
//             const integer dimmx)
//{
//    return ((y*dimmx)+x)*dimmz + z;
//};


#define STENCILZ(ptr,off,dzi,z,x,y,dimmz,dimmx) { \
   ((C0 * ( ptr[IDX(z  +off,x,y,dimmz,dimmx)] - ptr[IDX(z-1+off,x,y,dimmz,dimmx)]) +         \
     C1 * ( ptr[IDX(z+1+off,x,y,dimmz,dimmx)] - ptr[IDX(z-2+off,x,y,dimmz,dimmx)]) +         \
     C2 * ( ptr[IDX(z+2+off,x,y,dimmz,dimmx)] - ptr[IDX(z-3+off,x,y,dimmz,dimmx)]) +         \
     C3 * ( ptr[IDX(z+3+off,x,y,dimmz,dimmx)] - ptr[IDX(z-4+off,x,y,dimmz,dimmx)])) * dzi )  \
}


inline
real stencil_Z ( const integer off,
                 const real* restrict ptr,
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

#define STENCILX(ptr,off,dxi,z,x,y,dimmz,dimmx) { \
   ((C0 * ( ptr[IDX(z,x  +off,y,dimmz,dimmx)] - ptr[IDX(z,x-1+off,y,dimmz,dimmx)]) +         \
     C1 * ( ptr[IDX(z,x+1+off,y,dimmz,dimmx)] - ptr[IDX(z,x-2+off,y,dimmz,dimmx)]) +         \
     C2 * ( ptr[IDX(z,x+2+off,y,dimmz,dimmx)] - ptr[IDX(z,x-3+off,y,dimmz,dimmx)]) +         \
     C3 * ( ptr[IDX(z,x+3+off,y,dimmz,dimmx)] - ptr[IDX(z,x-4+off,y,dimmz,dimmx)])) * dxi )  \
}

inline
real stencil_X( const integer off,
                const real* restrict ptr,
                const real dxi,
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

#define STENCILY(ptr,off,dyi,z,x,y,dimmz,dimmx) { \
   ((C0 * ( ptr[IDX(z,x,y  +off,dimmz,dimmx)] - ptr[IDX(z,x,y-1+off,dimmz,dimmx)]) +         \
     C1 * ( ptr[IDX(z,x,y+1+off,dimmz,dimmx)] - ptr[IDX(z,x,y-2+off,dimmz,dimmx)]) +         \
     C2 * ( ptr[IDX(z,x,y+2+off,dimmz,dimmx)] - ptr[IDX(z,x,y-3+off,dimmz,dimmx)]) +         \
     C3 * ( ptr[IDX(z,x,y+3+off,dimmz,dimmx)] - ptr[IDX(z,x,y-4+off,dimmz,dimmx)])) * dyi )  \
}

inline
real stencil_Y( const integer off,
                const real* restrict ptr,
                const real dyi,
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

#define RHO_BL(ptr,z,x,y,dimmz,dimmx) { \
   (2.0f / (ptr[IDX(z,x,y,dimmz,dimmx)] + ptr[IDX(z+1,x,y,dimmz,dimmx)])) \
}

inline
real rho_BL ( const real* restrict rho,
              const integer z,
              const integer x,
              const integer y,
              const integer dimmz,
              const integer dimmx)
{
    return (2.0f / (rho[IDX(z,x,y,dimmz,dimmx)] + rho[IDX(z+1,x,y,dimmz,dimmx)]));
};

#define RHO_TR(ptr,z,x,y,dimmz,dimmx) { \
   (2.0f / (ptr[IDX(z,x,y,dimmz,dimmx)] + ptr[IDX(z,x+1,y,dimmz,dimmx)])) \
}

inline
real rho_TR ( const real* restrict rho,
              const integer z,
              const integer x,
              const integer y,
              const integer dimmz,
              const integer dimmx)
{
    return (2.0f / (rho[IDX(z,x,y,dimmz,dimmx)] + rho[IDX(z,x+1,y,dimmz,dimmx)]));
};

#define RHO_BR(ptr,z,x,y,dimmz,dimmx) { \
   ( 8.0f/ ( ptr[IDX(z  ,x  ,y  ,dimmz,dimmx)] + \
             ptr[IDX(z+1,x  ,y  ,dimmz,dimmx)] + \
             ptr[IDX(z  ,x+1,y  ,dimmz,dimmx)] + \
             ptr[IDX(z  ,x  ,y+1,dimmz,dimmx)] + \
             ptr[IDX(z  ,x+1,y+1,dimmz,dimmx)] + \
             ptr[IDX(z+1,x+1,y  ,dimmz,dimmx)] + \
             ptr[IDX(z+1,x  ,y+1,dimmz,dimmx)] + \
             ptr[IDX(z+1,x+1,y+1,dimmz,dimmx)])) \
}

inline
real rho_BR ( const real* restrict rho,
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

#define RHO_TL(ptr,z,x,y,dimmz,dimmx) { \
   (2.0f / (ptr[IDX(z,x,y,dimmz,dimmx)] + ptr[IDX(z,x,y+1,dimmz,dimmx)])) \
}

inline
real rho_TL ( const real* restrict rho,
              const integer z,
              const integer x,
              const integer y,
              const integer dimmz,
              const integer dimmx)
{
    return (2.0f / (rho[IDX(z,x,y,dimmz,dimmx)] + rho[IDX(z,x,y+1,dimmz,dimmx)]));
};

void compute_component_vcell_TL (      real* restrict vptr,
                                 const real* restrict szptr,
                                 const real* restrict sxptr,
                                 const real* restrict syptr,
                                 const real* restrict rho,
                                 const real           dt,
                                 const real           dzi,
                                 const real           dxi,
                                 const real           dyi,
                                 const integer        nz0,
                                 const integer        nzf,
                                 const integer        nx0,
                                 const integer        nxf,
                                 const integer        ny0,
                                 const integer        nyf,
                                 const offset_t       _SZ,
                                 const offset_t       _SX,
                                 const offset_t       _SY,
                                 const integer        dimmz,
                                 const integer        dimmx,
                                 const integer        dimmy,
                                 const phase_t        phase)
{
    const int SX = _SX;
    const int SY = _SY;
    const int SZ = _SZ;

    const integer size = dimmz * dimmx * dimmy;

    #pragma omp target device(openacc) copy_deps
    #pragma omp task in( [size]rho, [size]sxptr, [size]syptr, [size]szptr ) \
                     inout( [size]vptr )\
                     label(vcell_TL)
    #pragma acc kernels deviceptr(rho, sxptr, syptr, szptr, vptr)
    #pragma acc loop independent
    for(integer y=ny0; y < nyf; y++)
    {
        #pragma acc loop independent
        for(integer x=nx0; x < nxf; x++)
        {
            #pragma acc loop independent
            for(integer z=nz0; z < nzf; z++)
            {
                const real lrho = RHO_BL( rho,z,x,y,dimmz,dimmx );

                const real stx = STENCILX( sxptr,SX,dxi,z,x,y,dimmz,dimmx );
                const real sty = STENCILY( syptr,SY,dyi,z,x,y,dimmz,dimmx );
                const real stz = STENCILZ( szptr,SZ,dzi,z,x,y,dimmz,dimmx );

                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
};

void compute_component_vcell_TR (      real* restrict vptr,
                                 const real* restrict szptr,
                                 const real* restrict sxptr,
                                 const real* restrict syptr,
                                 const real* restrict rho,
                                 const real           dt,
                                 const real           dzi,
                                 const real           dxi,
                                 const real           dyi,
                                 const integer        nz0,
                                 const integer        nzf,
                                 const integer        nx0,
                                 const integer        nxf,
                                 const integer        ny0,
                                 const integer        nyf,
                                 const offset_t       _SZ,
                                 const offset_t       _SX,
                                 const offset_t       _SY,
                                 const integer        dimmz,
                                 const integer        dimmx,
                                 const integer        dimmy,
                                 const phase_t        phase)
{
    const int SX = _SX;
    const int SY = _SY;
    const int SZ = _SZ;

    const integer size = dimmz * dimmx * dimmy;

    #pragma omp target device(openacc) copy_deps
    #pragma omp task in( [size]rho, [size]sxptr, [size]syptr, [size]szptr ) \
                     inout( [size]vptr )\
                     label(vcell_TR)
    #pragma acc kernels deviceptr(rho, sxptr, syptr, szptr, vptr)
    #pragma acc loop independent
    for(integer y=ny0; y < nyf; y++)
    {
        #pragma acc loop independent
        for(integer x=nx0; x < nxf; x++)
        {
            #pragma acc loop independent
            for(integer z=nz0; z < nzf; z++)
            {
                const real lrho = RHO_TR( rho,z,x,y,dimmz,dimmx );

                const real stx = STENCILX( sxptr,SX,dxi,z,x,y,dimmz,dimmx );
                const real sty = STENCILY( syptr,SY,dyi,z,x,y,dimmz,dimmx );
                const real stz = STENCILZ( szptr,SZ,dzi,z,x,y,dimmz,dimmx );

                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
};

void compute_component_vcell_BR (      real* restrict vptr,
                                 const real* restrict szptr,
                                 const real* restrict sxptr,
                                 const real* restrict syptr,
                                 const real* restrict rho,
                                 const real           dt,
                                 const real           dzi,
                                 const real           dxi,
                                 const real           dyi,
                                 const integer        nz0,
                                 const integer        nzf,
                                 const integer        nx0,
                                 const integer        nxf,
                                 const integer        ny0,
                                 const integer        nyf,
                                 const offset_t       _SZ,
                                 const offset_t       _SX,
                                 const offset_t       _SY,
                                 const integer        dimmz,
                                 const integer        dimmx,
                                 const integer        dimmy,
                                 const phase_t        phase)
{
    const int SX = _SX;
    const int SY = _SY;
    const int SZ = _SZ;

    const integer size = dimmz * dimmx * dimmy;

    #pragma omp target device(openacc) copy_deps
    #pragma omp task in( [size]rho, [size]sxptr, [size]syptr, [size]szptr ) \
                     inout( [size]vptr )\
                     label(vcell_BR)
    #pragma acc kernels deviceptr(rho, sxptr, syptr, szptr, vptr)
    #pragma acc loop independent
    for(integer y=ny0; y < nyf; y++)
    {
        #pragma acc loop independent
        for(integer x=nx0; x < nxf; x++)
        {
            #pragma acc loop independent
            for(integer z=nz0; z < nzf; z++)
            {
                const real lrho = RHO_BR( rho,z,x,y,dimmz,dimmx );

                const real stx = STENCILX( sxptr,SX,dxi,z,x,y,dimmz,dimmx );
                const real sty = STENCILY( syptr,SY,dyi,z,x,y,dimmz,dimmx );
                const real stz = STENCILZ( szptr,SZ,dzi,z,x,y,dimmz,dimmx );

                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
};

void compute_component_vcell_BL (      real* restrict vptr,
                                 const real* restrict szptr,
                                 const real* restrict sxptr,
                                 const real* restrict syptr,
                                 const real* restrict rho,
                                 const real           dt,
                                 const real           dzi,
                                 const real           dxi,
                                 const real           dyi,
                                 const integer        nz0,
                                 const integer        nzf,
                                 const integer        nx0,
                                 const integer        nxf,
                                 const integer        ny0,
                                 const integer        nyf,
                                 const offset_t       _SZ,
                                 const offset_t       _SX,
                                 const offset_t       _SY,
                                 const integer        dimmz,
                                 const integer        dimmx,
                                 const integer        dimmy,
                                 const phase_t        phase)
{
    const int SX = _SX;
    const int SY = _SY;
    const int SZ = _SZ;

    const integer size = dimmz * dimmx * dimmy;

    #pragma omp target device(openacc) copy_deps
    #pragma omp task in( [size]rho, [size]sxptr, [size]syptr, [size]szptr ) \
                     inout( [size]vptr )\
                     label(vcell_BR)
    #pragma acc kernels deviceptr(rho, sxptr, syptr, szptr, vptr)
    #pragma acc loop independent
    for(integer y=ny0; y < nyf; y++)
    {
        #pragma acc loop independent
        for(integer x=nx0; x < nxf; x++)
        {
            #pragma acc loop independent
            for(integer z=nz0; z < nzf; z++)
            {
                const real lrho = RHO_BL( rho,z,x,y,dimmz,dimmx );

                const real stx = STENCILX( sxptr,SX,dxi,z,x,y,dimmz,dimmx );
                const real sty = STENCILY( syptr,SY,dyi,z,x,y,dimmz,dimmx );
                const real stz = STENCILZ( szptr,SZ,dzi,z,x,y,dimmz,dimmx );

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
                         const integer dimmy,
                         const phase_t phase)
{
#if defined(DEBUG)
    fprintf(stderr, "Integration limits of %s are (z "I"-"I",x "I"-"I",y "I"-"I")\n", __FUNCTION__, nz0,nzf,nx0,nxf,ny0,nyf);
#endif

#if defined(__INTEL_COMPILER)
    #pragma forceinline recursive
#endif
    {
        compute_component_vcell_TL (v.tl.w, s.bl.zz, s.tr.xz, s.tl.yz, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, back_offset, forw_offset, dimmz, dimmx, dimmy, phase);
        compute_component_vcell_TR (v.tr.w, s.br.zz, s.tl.xz, s.tr.yz, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, forw_offset, back_offset, dimmz, dimmx, dimmy, phase);
        compute_component_vcell_BL (v.bl.w, s.tl.zz, s.br.xz, s.bl.yz, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, back_offset, dimmz, dimmx, dimmy, phase);
        compute_component_vcell_BR (v.br.w, s.tr.zz, s.bl.xz, s.br.yz, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, forw_offset, forw_offset, dimmz, dimmx, dimmy, phase);
        compute_component_vcell_TL (v.tl.u, s.bl.xz, s.tr.xx, s.tl.xy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, back_offset, forw_offset, dimmz, dimmx, dimmy, phase);
        compute_component_vcell_TR (v.tr.u, s.br.xz, s.tl.xx, s.tr.xy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, forw_offset, back_offset, dimmz, dimmx, dimmy, phase);
        compute_component_vcell_BL (v.bl.u, s.tl.xz, s.br.xx, s.bl.xy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, back_offset, dimmz, dimmx, dimmy, phase);
        compute_component_vcell_BR (v.br.u, s.tr.xz, s.bl.xx, s.br.xy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, forw_offset, forw_offset, dimmz, dimmx, dimmy, phase);
        compute_component_vcell_TL (v.tl.v, s.bl.yz, s.tr.xy, s.tl.yy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, back_offset, forw_offset, dimmz, dimmx, dimmy, phase);
        compute_component_vcell_TR (v.tr.v, s.br.yz, s.tl.xy, s.tr.yy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, forw_offset, back_offset, dimmz, dimmx, dimmy, phase);
        compute_component_vcell_BL (v.bl.v, s.tl.yz, s.br.xy, s.bl.yy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, back_offset, dimmz, dimmx, dimmy, phase);
        compute_component_vcell_BR (v.br.v, s.tr.yz, s.bl.xy, s.br.yy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, forw_offset, forw_offset, dimmz, dimmx, dimmy, phase);
    }
};





/* ------------------------------------------------------------------------------ */
/*                                                                                */
/*                               CALCULO DE TENSIONES                             */
/*                                                                                */
/* ------------------------------------------------------------------------------ */

#define STRESS_UPDATE(sptr,c1,c2,c3,c4,c5,c6,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx) { \
   {                                           \
      real accum  = dt * c1 * u_x;             \
           accum += dt * c2 * v_y;             \
           accum += dt * c3 * w_z;             \
           accum += dt * c4 * (w_y + v_z);     \
           accum += dt * c5 * (w_x + u_z);     \
           accum += dt * c6 * (v_x + u_y);     \
      sptr[IDX(z,x,y,dimmz,dimmx)] += accum;   \
   }                                           \
}

inline
void stress_update(real* restrict sptr,
                   const real     c1,
                   const real     c2,
                   const real     c3,
                   const real     c4,
                   const real     c5,
                   const real     c6,
                   const integer  z,
                   const integer  x,
                   const integer  y,
                   const real     dt,
                   const real     u_x,
                   const real     u_y,
                   const real     u_z,
                   const real     v_x,
                   const real     v_y,
                   const real     v_z,
                   const real     w_x,
                   const real     w_y,
                   const real     w_z,
                   const integer  dimmz,
                   const integer  dimmx)
{
    real accum  = dt * c1 * u_x;
         accum += dt * c2 * v_y;
         accum += dt * c3 * w_z;
         accum += dt * c4 * (w_y + v_z);
         accum += dt * c5 * (w_x + u_z);
         accum += dt * c6 * (v_x + u_y);
    sptr[IDX(z,x,y,dimmz,dimmx)] += accum;
};

#define CELL_COEFF_BR(ptr,z,x,y,dimmz,dimmx) {            \
   ( 1.0f / ( 2.5f  *(ptr[IDX(z  ,x  ,y,dimmz,dimmx)] +  \
                      ptr[IDX(z  ,x+1,y,dimmz,dimmx)] +  \
                      ptr[IDX(z+1,x  ,y,dimmz,dimmx)] +  \
                      ptr[IDX(z+1,x+1,y,dimmz,dimmx)])) )\
}

inline
real cell_coeff_BR ( const real* restrict ptr,
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

#define CELL_COEFF_TL(ptr,z,x,y,dimmz,dimmx) { \
    ( 1.0f / (ptr[IDX(z,x,y,dimmz,dimmx)]) )   \
}

inline
real cell_coeff_TL ( const real* restrict ptr,
                     const integer z,
                     const integer x,
                     const integer y,
                     const integer dimmz,
                     const integer dimmx)
{
    return ( 1.0f / (ptr[IDX(z,x,y,dimmz,dimmx)]));
};

#define CELL_COEFF_BL(ptr,z,x,y,dimmz,dimmx) {            \
   ( 1.0f / ( 2.5f  *(ptr[IDX(z  ,x,y  ,dimmz,dimmx)] +  \
                      ptr[IDX(z  ,x,y+1,dimmz,dimmx)] +  \
                      ptr[IDX(z+1,x,y  ,dimmz,dimmx)] +  \
                      ptr[IDX(z+1,x,y+1,dimmz,dimmx)])) )\
}

inline
real cell_coeff_BL ( const real* restrict ptr,
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

#define CELL_COEFF_TR(ptr,z,x,y,dimmz,dimmx) {            \
   ( 1.0f / ( 2.5f  *(ptr[IDX(z,x  ,y  ,dimmz,dimmx)] +  \
                      ptr[IDX(z,x+1,y  ,dimmz,dimmx)] +  \
                      ptr[IDX(z,x  ,y+1,dimmz,dimmx)] +  \
                      ptr[IDX(z,x+1,y+1,dimmz,dimmx)])) )\
}

inline
real cell_coeff_TR ( const real* restrict ptr,
                     const integer z,
                     const integer x,
                     const integer y,
                     const integer dimmz,
                     const integer dimmx)
{
    return ( 1.0f / ( 2.5f *(ptr[IDX(z,x  , y  ,dimmz,dimmx)] +
                             ptr[IDX(z,x+1, y  ,dimmz,dimmx)] +
                             ptr[IDX(z,x  , y+1,dimmz,dimmx)] +
                             ptr[IDX(z,x+1, y+1,dimmz,dimmx)])));
};

#define CELL_COEFF_ARTM_BR(ptr,z,x,y,dimmz,dimmx) {     \
    ((1.0f / ptr[IDX(z  ,x  ,y,dimmz,dimmx )]  +        \
      1.0f / ptr[IDX(z  ,x+1,y,dimmz,dimmx )]  +        \
      1.0f / ptr[IDX(z+1,x  ,y,dimmz,dimmx )]  +        \
      1.0f / ptr[IDX(z+1,x+1,y,dimmz,dimmx )]) * 0.25f) \
}

inline
real cell_coeff_ARTM_BR( const real* restrict ptr,
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

#define CELL_COEFF_ARTM_TL(ptr,z,x,y,dimmz,dimmx) { \
    ( 1.0f / (ptr[IDX(z,x,y,dimmz,dimmx)]) )   \
}

inline
real cell_coeff_ARTM_TL( const real* restrict ptr,
                         const integer z,
                         const integer x,
                         const integer y,
                         const integer dimmz,
                         const integer dimmx)
{
    return (1.0f / ptr[IDX(z,x,y,dimmz,dimmx)]);
};

#define CELL_COEFF_ARTM_BL(ptr,z,x,y,dimmz,dimmx) {     \
    ((1.0f / ptr[IDX(z  ,x,y  ,dimmz,dimmx )]  +        \
      1.0f / ptr[IDX(z  ,x,y+1,dimmz,dimmx )]  +        \
      1.0f / ptr[IDX(z+1,x,y  ,dimmz,dimmx )]  +        \
      1.0f / ptr[IDX(z+1,x,y+1,dimmz,dimmx )]) * 0.25f) \
}

inline
real cell_coeff_ARTM_BL( const real* restrict ptr,
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

#define CELL_COEFF_ARTM_TR(ptr,z,x,y,dimmz,dimmx) {     \
    ((1.0f / ptr[IDX(z,x  ,y  ,dimmz,dimmx )]  +        \
      1.0f / ptr[IDX(z,x+1,y  ,dimmz,dimmx )]  +        \
      1.0f / ptr[IDX(z,x  ,y+1,dimmz,dimmx )]  +        \
      1.0f / ptr[IDX(z,x+1,y+1,dimmz,dimmx )]) * 0.25f) \
}

inline
real cell_coeff_ARTM_TR( const real* restrict ptr,
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
                                 const integer  dimmx,
                                 const integer  dimmy,
                                 const phase_t phase)
{
    real* restrict sxxptr __attribute__ ((aligned (64))) = s.tr.xx;
    real* restrict syyptr __attribute__ ((aligned (64))) = s.tr.yy;
    real* restrict szzptr __attribute__ ((aligned (64))) = s.tr.zz;
    real* restrict syzptr __attribute__ ((aligned (64))) = s.tr.yz;
    real* restrict sxzptr __attribute__ ((aligned (64))) = s.tr.xz;
    real* restrict sxyptr __attribute__ ((aligned (64))) = s.tr.xy;

    const real* restrict vxu    __attribute__ ((aligned (64))) = vnode_x.u;
    const real* restrict vxv    __attribute__ ((aligned (64))) = vnode_x.v;
    const real* restrict vxw    __attribute__ ((aligned (64))) = vnode_x.w;
    const real* restrict vyu    __attribute__ ((aligned (64))) = vnode_y.u;
    const real* restrict vyv    __attribute__ ((aligned (64))) = vnode_y.v;
    const real* restrict vyw    __attribute__ ((aligned (64))) = vnode_y.w;
    const real* restrict vzu    __attribute__ ((aligned (64))) = vnode_z.u;
    const real* restrict vzv    __attribute__ ((aligned (64))) = vnode_z.v;
    const real* restrict vzw    __attribute__ ((aligned (64))) = vnode_z.w;

    const real* restrict cc11 = coeffs.c11;
    const real* restrict cc12 = coeffs.c12;
    const real* restrict cc13 = coeffs.c13;
    const real* restrict cc14 = coeffs.c14;
    const real* restrict cc15 = coeffs.c15;
    const real* restrict cc16 = coeffs.c16;
    const real* restrict cc22 = coeffs.c22;
    const real* restrict cc23 = coeffs.c23;
    const real* restrict cc24 = coeffs.c24;
    const real* restrict cc25 = coeffs.c25;
    const real* restrict cc26 = coeffs.c26;
    const real* restrict cc33 = coeffs.c33;
    const real* restrict cc34 = coeffs.c34;
    const real* restrict cc35 = coeffs.c35;
    const real* restrict cc36 = coeffs.c36;
    const real* restrict cc44 = coeffs.c44;
    const real* restrict cc45 = coeffs.c45;
    const real* restrict cc46 = coeffs.c46;
    const real* restrict cc55 = coeffs.c55;
    const real* restrict cc56 = coeffs.c56;
    const real* restrict cc66 = coeffs.c66;

    const int SX = _SX;
    const int SY = _SY;
    const int SZ = _SZ;

    const integer size = dimmz * dimmx * dimmy;

    #pragma omp target device(openacc) copy_deps
    #pragma omp task in( [size]cc11, [size]cc12, [size]cc13, [size]cc14, [size]cc15, [size]cc16 ) \
                     in(             [size]cc22, [size]cc23, [size]cc24, [size]cc25, [size]cc26 ) \
                     in(                         [size]cc33, [size]cc34, [size]cc35, [size]cc36 ) \
                     in(                                     [size]cc44, [size]cc45, [size]cc46 ) \
                     in(                                                 [size]cc55, [size]cc56 ) \
                     in(                                                             [size]cc66 ) \
                     in( [size]vxu, [size]vxv, [size]vxw, [size]vyu, [size]vyv, [size]vyw, [size]vzu, [size]vzv, [size]vzw ) \
                     inout( [size]sxxptr, [size]syyptr, [size]szzptr, [size]syzptr, [size]sxzptr, [size]sxyptr ) \
                     label(scell_TR)
    #pragma acc kernels deviceptr(cc11, cc12, cc13, cc14, cc15, cc16) \
                        deviceptr(      cc22, cc23, cc24, cc25, cc26) \
                        deviceptr(            cc33, cc34, cc35, cc36) \
                        deviceptr(                  cc44, cc45, cc46) \
                        deviceptr(                        cc55, cc56) \
                        deviceptr(                              cc66) \
                        deviceptr(vxu, vxv, vxw) \
                        deviceptr(vyu, vyv, vyw) \
                        deviceptr(vzu, vzv, vzw) \
                        deviceptr(sxxptr, syyptr, szzptr, syzptr, sxzptr, sxyptr)
    #pragma acc loop independent
    for (integer y = ny0; y < nyf; y++)
    {
        #pragma acc loop independent
        for (integer x = nx0; x < nxf; x++)
        {
            #pragma acc loop independent
            for (integer z = nz0; z < nzf; z++ )
            {
                const real c11 = CELL_COEFF_TR      (cc11, z, x, y, dimmz, dimmx);
                const real c12 = CELL_COEFF_TR      (cc12, z, x, y, dimmz, dimmx);
                const real c13 = CELL_COEFF_TR      (cc13, z, x, y, dimmz, dimmx);
                const real c14 = CELL_COEFF_ARTM_TR (cc14, z, x, y, dimmz, dimmx);
                const real c15 = CELL_COEFF_ARTM_TR (cc15, z, x, y, dimmz, dimmx);
                const real c16 = CELL_COEFF_ARTM_TR (cc16, z, x, y, dimmz, dimmx);
                const real c22 = CELL_COEFF_TR      (cc22, z, x, y, dimmz, dimmx);
                const real c23 = CELL_COEFF_TR      (cc23, z, x, y, dimmz, dimmx);
                const real c24 = CELL_COEFF_ARTM_TR (cc24, z, x, y, dimmz, dimmx);
                const real c25 = CELL_COEFF_ARTM_TR (cc25, z, x, y, dimmz, dimmx);
                const real c26 = CELL_COEFF_ARTM_TR (cc26, z, x, y, dimmz, dimmx);
                const real c33 = CELL_COEFF_TR      (cc33, z, x, y, dimmz, dimmx);
                const real c34 = CELL_COEFF_ARTM_TR (cc34, z, x, y, dimmz, dimmx);
                const real c35 = CELL_COEFF_ARTM_TR (cc35, z, x, y, dimmz, dimmx);
                const real c36 = CELL_COEFF_ARTM_TR (cc36, z, x, y, dimmz, dimmx);
                const real c44 = CELL_COEFF_TR      (cc44, z, x, y, dimmz, dimmx);
                const real c45 = CELL_COEFF_ARTM_TR (cc45, z, x, y, dimmz, dimmx);
                const real c46 = CELL_COEFF_ARTM_TR (cc46, z, x, y, dimmz, dimmx);
                const real c55 = CELL_COEFF_TR      (cc55, z, x, y, dimmz, dimmx);
                const real c56 = CELL_COEFF_ARTM_TR (cc56, z, x, y, dimmz, dimmx);
                const real c66 = CELL_COEFF_TR      (cc66, z, x, y, dimmz, dimmx);

                const real u_x = STENCILX (vxu, SX, dxi, z, x, y, dimmz, dimmx);
                const real v_x = STENCILX (vxv, SX, dxi, z, x, y, dimmz, dimmx);
                const real w_x = STENCILX (vxw, SX, dxi, z, x, y, dimmz, dimmx);

                const real u_y = STENCILY (vyu, SY, dyi, z, x, y, dimmz, dimmx);
                const real v_y = STENCILY (vyv, SY, dyi, z, x, y, dimmz, dimmx);
                const real w_y = STENCILY (vyw, SY, dyi, z, x, y, dimmz, dimmx);

                const real u_z = STENCILZ (vzu, SZ, dzi, z, x, y, dimmz, dimmx);
                const real v_z = STENCILZ (vzv, SZ, dzi, z, x, y, dimmz, dimmx);
                const real w_z = STENCILZ (vzw, SZ, dzi, z, x, y, dimmz, dimmx);

                STRESS_UPDATE (sxxptr,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                STRESS_UPDATE (syyptr,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                STRESS_UPDATE (szzptr,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                STRESS_UPDATE (syzptr,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                STRESS_UPDATE (sxzptr,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                STRESS_UPDATE (sxyptr,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            }
        }
    }
};

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
                                 const integer  dimmx,
                                 const integer  dimmy,
                                 const phase_t phase)
{
    real* restrict sxxptr __attribute__ ((aligned (64))) = s.tl.xx;
    real* restrict syyptr __attribute__ ((aligned (64))) = s.tl.yy;
    real* restrict szzptr __attribute__ ((aligned (64))) = s.tl.zz;
    real* restrict syzptr __attribute__ ((aligned (64))) = s.tl.yz;
    real* restrict sxzptr __attribute__ ((aligned (64))) = s.tl.xz;
    real* restrict sxyptr __attribute__ ((aligned (64))) = s.tl.xy;

    const real* restrict vxu    __attribute__ ((aligned (64))) = vnode_x.u;
    const real* restrict vxv    __attribute__ ((aligned (64))) = vnode_x.v;
    const real* restrict vxw    __attribute__ ((aligned (64))) = vnode_x.w;
    const real* restrict vyu    __attribute__ ((aligned (64))) = vnode_y.u;
    const real* restrict vyv    __attribute__ ((aligned (64))) = vnode_y.v;
    const real* restrict vyw    __attribute__ ((aligned (64))) = vnode_y.w;
    const real* restrict vzu    __attribute__ ((aligned (64))) = vnode_z.u;
    const real* restrict vzv    __attribute__ ((aligned (64))) = vnode_z.v;
    const real* restrict vzw    __attribute__ ((aligned (64))) = vnode_z.w;

    const real* restrict cc11 = coeffs.c11;
    const real* restrict cc12 = coeffs.c12;
    const real* restrict cc13 = coeffs.c13;
    const real* restrict cc14 = coeffs.c14;
    const real* restrict cc15 = coeffs.c15;
    const real* restrict cc16 = coeffs.c16;
    const real* restrict cc22 = coeffs.c22;
    const real* restrict cc23 = coeffs.c23;
    const real* restrict cc24 = coeffs.c24;
    const real* restrict cc25 = coeffs.c25;
    const real* restrict cc26 = coeffs.c26;
    const real* restrict cc33 = coeffs.c33;
    const real* restrict cc34 = coeffs.c34;
    const real* restrict cc35 = coeffs.c35;
    const real* restrict cc36 = coeffs.c36;
    const real* restrict cc44 = coeffs.c44;
    const real* restrict cc45 = coeffs.c45;
    const real* restrict cc46 = coeffs.c46;
    const real* restrict cc55 = coeffs.c55;
    const real* restrict cc56 = coeffs.c56;
    const real* restrict cc66 = coeffs.c66;

    const int SX = _SX;
    const int SY = _SY;
    const int SZ = _SZ;

    const integer size = dimmz * dimmx * dimmy;

    #pragma omp target device(openacc) copy_deps
    #pragma omp task in( [size]cc11, [size]cc12, [size]cc13, [size]cc14, [size]cc15, [size]cc16 ) \
                     in(             [size]cc22, [size]cc23, [size]cc24, [size]cc25, [size]cc26 ) \
                     in(                         [size]cc33, [size]cc34, [size]cc35, [size]cc36 ) \
                     in(                                     [size]cc44, [size]cc45, [size]cc46 ) \
                     in(                                                 [size]cc55, [size]cc56 ) \
                     in(                                                             [size]cc66 ) \
                     in( [size]vxu, [size]vxv, [size]vxw, [size]vyu, [size]vyv, [size]vyw, [size]vzu, [size]vzv, [size]vzw ) \
                     inout( [size]sxxptr, [size]syyptr, [size]szzptr, [size]syzptr, [size]sxzptr, [size]sxyptr ) \
                     label(scell_TL)
    #pragma acc kernels deviceptr(cc11, cc12, cc13, cc14, cc15, cc16) \
                        deviceptr(      cc22, cc23, cc24, cc25, cc26) \
                        deviceptr(            cc33, cc34, cc35, cc36) \
                        deviceptr(                  cc44, cc45, cc46) \
                        deviceptr(                        cc55, cc56) \
                        deviceptr(                              cc66) \
                        deviceptr(vxu, vxv, vxw) \
                        deviceptr(vyu, vyv, vyw) \
                        deviceptr(vzu, vzv, vzw) \
                        deviceptr(sxxptr, syyptr, szzptr, syzptr, sxzptr, sxyptr)
    #pragma acc loop independent
    for (integer y = ny0; y < nyf; y++)
    {
        #pragma acc loop independent
        for (integer x = nx0; x < nxf; x++)
        {
            #pragma acc loop independent
            for (integer z = nz0; z < nzf; z++ )
            {
                const real c11 = CELL_COEFF_TL      (cc11, z, x, y, dimmz, dimmx);
                const real c12 = CELL_COEFF_TL      (cc12, z, x, y, dimmz, dimmx);
                const real c13 = CELL_COEFF_TL      (cc13, z, x, y, dimmz, dimmx);
                const real c14 = CELL_COEFF_ARTM_TL (cc14, z, x, y, dimmz, dimmx);
                const real c15 = CELL_COEFF_ARTM_TL (cc15, z, x, y, dimmz, dimmx);
                const real c16 = CELL_COEFF_ARTM_TL (cc16, z, x, y, dimmz, dimmx);
                const real c22 = CELL_COEFF_TL      (cc22, z, x, y, dimmz, dimmx);
                const real c23 = CELL_COEFF_TL      (cc23, z, x, y, dimmz, dimmx);
                const real c24 = CELL_COEFF_ARTM_TL (cc24, z, x, y, dimmz, dimmx);
                const real c25 = CELL_COEFF_ARTM_TL (cc25, z, x, y, dimmz, dimmx);
                const real c26 = CELL_COEFF_ARTM_TL (cc26, z, x, y, dimmz, dimmx);
                const real c33 = CELL_COEFF_TL      (cc33, z, x, y, dimmz, dimmx);
                const real c34 = CELL_COEFF_ARTM_TL (cc34, z, x, y, dimmz, dimmx);
                const real c35 = CELL_COEFF_ARTM_TL (cc35, z, x, y, dimmz, dimmx);
                const real c36 = CELL_COEFF_ARTM_TL (cc36, z, x, y, dimmz, dimmx);
                const real c44 = CELL_COEFF_TL      (cc44, z, x, y, dimmz, dimmx);
                const real c45 = CELL_COEFF_ARTM_TL (cc45, z, x, y, dimmz, dimmx);
                const real c46 = CELL_COEFF_ARTM_TL (cc46, z, x, y, dimmz, dimmx);
                const real c55 = CELL_COEFF_TL      (cc55, z, x, y, dimmz, dimmx);
                const real c56 = CELL_COEFF_ARTM_TL (cc56, z, x, y, dimmz, dimmx);
                const real c66 = CELL_COEFF_TL      (cc66, z, x, y, dimmz, dimmx);

                const real u_x = STENCILX (vxu, SX, dxi, z, x, y, dimmz, dimmx);
                const real v_x = STENCILX (vxv, SX, dxi, z, x, y, dimmz, dimmx);
                const real w_x = STENCILX (vxw, SX, dxi, z, x, y, dimmz, dimmx);

                const real u_y = STENCILY (vyu, SY, dyi, z, x, y, dimmz, dimmx);
                const real v_y = STENCILY (vyv, SY, dyi, z, x, y, dimmz, dimmx);
                const real w_y = STENCILY (vyw, SY, dyi, z, x, y, dimmz, dimmx);

                const real u_z = STENCILZ (vzu, SZ, dzi, z, x, y, dimmz, dimmx);
                const real v_z = STENCILZ (vzv, SZ, dzi, z, x, y, dimmz, dimmx);
                const real w_z = STENCILZ (vzw, SZ, dzi, z, x, y, dimmz, dimmx);

                STRESS_UPDATE (sxxptr,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                STRESS_UPDATE (syyptr,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                STRESS_UPDATE (szzptr,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                STRESS_UPDATE (syzptr,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                STRESS_UPDATE (sxzptr,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                STRESS_UPDATE (sxyptr,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            }
        }
    }
};


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
                                 const integer  dimmx,
                                 const integer  dimmy,
                                 const phase_t phase)
{
    real* restrict sxxptr __attribute__ ((aligned (64))) = s.br.xx;
    real* restrict syyptr __attribute__ ((aligned (64))) = s.br.yy;
    real* restrict szzptr __attribute__ ((aligned (64))) = s.br.zz;
    real* restrict syzptr __attribute__ ((aligned (64))) = s.br.yz;
    real* restrict sxzptr __attribute__ ((aligned (64))) = s.br.xz;
    real* restrict sxyptr __attribute__ ((aligned (64))) = s.br.xy;

    const real* restrict vxu    __attribute__ ((aligned (64))) = vnode_x.u;
    const real* restrict vxv    __attribute__ ((aligned (64))) = vnode_x.v;
    const real* restrict vxw    __attribute__ ((aligned (64))) = vnode_x.w;
    const real* restrict vyu    __attribute__ ((aligned (64))) = vnode_y.u;
    const real* restrict vyv    __attribute__ ((aligned (64))) = vnode_y.v;
    const real* restrict vyw    __attribute__ ((aligned (64))) = vnode_y.w;
    const real* restrict vzu    __attribute__ ((aligned (64))) = vnode_z.u;
    const real* restrict vzv    __attribute__ ((aligned (64))) = vnode_z.v;
    const real* restrict vzw    __attribute__ ((aligned (64))) = vnode_z.w;

    const real* restrict cc11 = coeffs.c11;
    const real* restrict cc12 = coeffs.c12;
    const real* restrict cc13 = coeffs.c13;
    const real* restrict cc14 = coeffs.c14;
    const real* restrict cc15 = coeffs.c15;
    const real* restrict cc16 = coeffs.c16;
    const real* restrict cc22 = coeffs.c22;
    const real* restrict cc23 = coeffs.c23;
    const real* restrict cc24 = coeffs.c24;
    const real* restrict cc25 = coeffs.c25;
    const real* restrict cc26 = coeffs.c26;
    const real* restrict cc33 = coeffs.c33;
    const real* restrict cc34 = coeffs.c34;
    const real* restrict cc35 = coeffs.c35;
    const real* restrict cc36 = coeffs.c36;
    const real* restrict cc44 = coeffs.c44;
    const real* restrict cc45 = coeffs.c45;
    const real* restrict cc46 = coeffs.c46;
    const real* restrict cc55 = coeffs.c55;
    const real* restrict cc56 = coeffs.c56;
    const real* restrict cc66 = coeffs.c66;

    const int SX = _SX;
    const int SY = _SY;
    const int SZ = _SZ;

    const integer size = dimmz * dimmx * dimmy;

    #pragma omp target device(openacc) copy_deps
    #pragma omp task in( [size]cc11, [size]cc12, [size]cc13, [size]cc14, [size]cc15, [size]cc16) \
                     in(             [size]cc22, [size]cc23, [size]cc24, [size]cc25, [size]cc26) \
                     in(                         [size]cc33, [size]cc34, [size]cc35, [size]cc36) \
                     in(                                     [size]cc44, [size]cc45, [size]cc46) \
                     in(                                                 [size]cc55, [size]cc56) \
                     in(                                                             [size]cc66) \
                     in( [size]vxu, [size]vxv, [size]vxw, [size]vyu, [size]vyv, [size]vyw, [size]vzu, [size]vzv, [size]vzw) \
                     inout( [size]sxxptr, [size]syyptr, [size]szzptr, [size]syzptr, [size]sxzptr, [size]sxyptr) \
                     label(scell_BR)
    #pragma acc kernels deviceptr(cc11, cc12, cc13, cc14, cc15, cc16) \
                        deviceptr(      cc22, cc23, cc24, cc25, cc26) \
                        deviceptr(            cc33, cc34, cc35, cc36) \
                        deviceptr(                  cc44, cc45, cc46) \
                        deviceptr(                        cc55, cc56) \
                        deviceptr(                              cc66) \
                        deviceptr(vxu, vxv, vxw) \
                        deviceptr(vyu, vyv, vyw) \
                        deviceptr(vzu, vzv, vzw) \
                        deviceptr(sxxptr, syyptr, szzptr, syzptr, sxzptr, sxyptr)
    #pragma acc loop independent
    for (integer y = ny0; y < nyf; y++)
    {
        #pragma acc loop independent
        for (integer x = nx0; x < nxf; x++)
        {
            #pragma acc loop independent
            for (integer z = nz0; z < nzf; z++ )
            {
                const real c11 = CELL_COEFF_BR      (cc11, z, x, y, dimmz, dimmx);
                const real c12 = CELL_COEFF_BR      (cc12, z, x, y, dimmz, dimmx);
                const real c13 = CELL_COEFF_BR      (cc13, z, x, y, dimmz, dimmx);
                const real c22 = CELL_COEFF_BR      (cc22, z, x, y, dimmz, dimmx);
                const real c23 = CELL_COEFF_BR      (cc23, z, x, y, dimmz, dimmx);
                const real c33 = CELL_COEFF_BR      (cc33, z, x, y, dimmz, dimmx);
                const real c44 = CELL_COEFF_BR      (cc44, z, x, y, dimmz, dimmx);
                const real c55 = CELL_COEFF_BR      (cc55, z, x, y, dimmz, dimmx);
                const real c66 = CELL_COEFF_BR      (cc66, z, x, y, dimmz, dimmx);

                const real c14 = CELL_COEFF_ARTM_BR (cc14, z, x, y, dimmz, dimmx);
                const real c15 = CELL_COEFF_ARTM_BR (cc15, z, x, y, dimmz, dimmx);
                const real c16 = CELL_COEFF_ARTM_BR (cc16, z, x, y, dimmz, dimmx);
                const real c24 = CELL_COEFF_ARTM_BR (cc24, z, x, y, dimmz, dimmx);
                const real c25 = CELL_COEFF_ARTM_BR (cc25, z, x, y, dimmz, dimmx);
                const real c26 = CELL_COEFF_ARTM_BR (cc26, z, x, y, dimmz, dimmx);
                const real c34 = CELL_COEFF_ARTM_BR (cc34, z, x, y, dimmz, dimmx);
                const real c35 = CELL_COEFF_ARTM_BR (cc35, z, x, y, dimmz, dimmx);
                const real c36 = CELL_COEFF_ARTM_BR (cc36, z, x, y, dimmz, dimmx);
                const real c45 = CELL_COEFF_ARTM_BR (cc45, z, x, y, dimmz, dimmx);
                const real c46 = CELL_COEFF_ARTM_BR (cc46, z, x, y, dimmz, dimmx);
                const real c56 = CELL_COEFF_ARTM_BR (cc56, z, x, y, dimmz, dimmx);

                const real u_x = STENCILX (vxu, SX, dxi, z, x, y, dimmz, dimmx);
                const real v_x = STENCILX (vxv, SX, dxi, z, x, y, dimmz, dimmx);
                const real w_x = STENCILX (vxw, SX, dxi, z, x, y, dimmz, dimmx);

                const real u_y = STENCILY (vyu, SY, dyi, z, x, y, dimmz, dimmx);
                const real v_y = STENCILY (vyv, SY, dyi, z, x, y, dimmz, dimmx);
                const real w_y = STENCILY (vyw, SY, dyi, z, x, y, dimmz, dimmx);

                const real u_z = STENCILZ (vzu, SZ, dzi, z, x, y, dimmz, dimmx);
                const real v_z = STENCILZ (vzv, SZ, dzi, z, x, y, dimmz, dimmx);
                const real w_z = STENCILZ (vzw, SZ, dzi, z, x, y, dimmz, dimmx);

                STRESS_UPDATE (sxxptr,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                STRESS_UPDATE (syyptr,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                STRESS_UPDATE (szzptr,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                STRESS_UPDATE (syzptr,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                STRESS_UPDATE (sxzptr,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                STRESS_UPDATE (sxyptr,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            }
        }
    }
};

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
                                 const integer  dimmx,
                                 const integer   dimmy,
                                 const phase_t phase)
{
    real* restrict sxxptr __attribute__ ((aligned (64))) = s.br.xx;
    real* restrict syyptr __attribute__ ((aligned (64))) = s.br.yy;
    real* restrict szzptr __attribute__ ((aligned (64))) = s.br.zz;
    real* restrict syzptr __attribute__ ((aligned (64))) = s.br.yz;
    real* restrict sxzptr __attribute__ ((aligned (64))) = s.br.xz;
    real* restrict sxyptr __attribute__ ((aligned (64))) = s.br.xy;

    const real* restrict vxu    __attribute__ ((aligned (64))) = vnode_x.u;
    const real* restrict vxv    __attribute__ ((aligned (64))) = vnode_x.v;
    const real* restrict vxw    __attribute__ ((aligned (64))) = vnode_x.w;
    const real* restrict vyu    __attribute__ ((aligned (64))) = vnode_y.u;
    const real* restrict vyv    __attribute__ ((aligned (64))) = vnode_y.v;
    const real* restrict vyw    __attribute__ ((aligned (64))) = vnode_y.w;
    const real* restrict vzu    __attribute__ ((aligned (64))) = vnode_z.u;
    const real* restrict vzv    __attribute__ ((aligned (64))) = vnode_z.v;
    const real* restrict vzw    __attribute__ ((aligned (64))) = vnode_z.w;

    const real* restrict cc11 = coeffs.c11;
    const real* restrict cc12 = coeffs.c12;
    const real* restrict cc13 = coeffs.c13;
    const real* restrict cc14 = coeffs.c14;
    const real* restrict cc15 = coeffs.c15;
    const real* restrict cc16 = coeffs.c16;
    const real* restrict cc22 = coeffs.c22;
    const real* restrict cc23 = coeffs.c23;
    const real* restrict cc24 = coeffs.c24;
    const real* restrict cc25 = coeffs.c25;
    const real* restrict cc26 = coeffs.c26;
    const real* restrict cc33 = coeffs.c33;
    const real* restrict cc34 = coeffs.c34;
    const real* restrict cc35 = coeffs.c35;
    const real* restrict cc36 = coeffs.c36;
    const real* restrict cc44 = coeffs.c44;
    const real* restrict cc45 = coeffs.c45;
    const real* restrict cc46 = coeffs.c46;
    const real* restrict cc55 = coeffs.c55;
    const real* restrict cc56 = coeffs.c56;
    const real* restrict cc66 = coeffs.c66;

    const int SX = _SX;
    const int SY = _SY;
    const int SZ = _SZ;

    const integer size = dimmz * dimmx * dimmy;

    #pragma omp target device(openacc) copy_deps
    #pragma omp task in( [size]cc11, [size]cc12, [size]cc13, [size]cc14, [size]cc15, [size]cc16) \
                     in(             [size]cc22, [size]cc23, [size]cc24, [size]cc25, [size]cc26) \
                     in(                         [size]cc33, [size]cc34, [size]cc35, [size]cc36) \
                     in(                                     [size]cc44, [size]cc45, [size]cc46) \
                     in(                                                 [size]cc55, [size]cc56) \
                     in(                                                             [size]cc66) \
                     in( [size]vxu, [size]vxv, [size]vxw, [size]vyu, [size]vyv, [size]vyw, [size]vzu, [size]vzv, [size]vzw) \
                     inout( [size]sxxptr, [size]syyptr, [size]szzptr, [size]syzptr, [size]sxzptr, [size]sxyptr) \
                     label(scell_BL)
    #pragma acc kernels deviceptr(cc11, cc12, cc13, cc14, cc15, cc16) \
                        deviceptr(      cc22, cc23, cc24, cc25, cc26) \
                        deviceptr(            cc33, cc34, cc35, cc36) \
                        deviceptr(                  cc44, cc45, cc46) \
                        deviceptr(                        cc55, cc56) \
                        deviceptr(                              cc66) \
                        deviceptr(vxu, vxv, vxw) \
                        deviceptr(vyu, vyv, vyw) \
                        deviceptr(vzu, vzv, vzw) \
                        deviceptr(sxxptr, syyptr, szzptr, syzptr, sxzptr, sxyptr)
    #pragma acc loop independent
    for (integer y = ny0; y < nyf; y++)
    {
        #pragma acc loop independent
        for (integer x = nx0; x < nxf; x++)
        {
            #pragma acc loop independent
            for (integer z = nz0; z < nzf; z++ )
            {
                const real c11 = CELL_COEFF_BL      (cc11, z, x, y, dimmz, dimmx);
                const real c12 = CELL_COEFF_BL      (cc12, z, x, y, dimmz, dimmx);
                const real c13 = CELL_COEFF_BL      (cc13, z, x, y, dimmz, dimmx);
                const real c14 = CELL_COEFF_ARTM_BL (cc14, z, x, y, dimmz, dimmx);
                const real c15 = CELL_COEFF_ARTM_BL (cc15, z, x, y, dimmz, dimmx);
                const real c16 = CELL_COEFF_ARTM_BL (cc16, z, x, y, dimmz, dimmx);
                const real c22 = CELL_COEFF_BL      (cc22, z, x, y, dimmz, dimmx);
                const real c23 = CELL_COEFF_BL      (cc23, z, x, y, dimmz, dimmx);
                const real c24 = CELL_COEFF_ARTM_BL (cc24, z, x, y, dimmz, dimmx);
                const real c25 = CELL_COEFF_ARTM_BL (cc25, z, x, y, dimmz, dimmx);
                const real c26 = CELL_COEFF_ARTM_BL (cc26, z, x, y, dimmz, dimmx);
                const real c33 = CELL_COEFF_BL      (cc33, z, x, y, dimmz, dimmx);
                const real c34 = CELL_COEFF_ARTM_BL (cc34, z, x, y, dimmz, dimmx);
                const real c35 = CELL_COEFF_ARTM_BL (cc35, z, x, y, dimmz, dimmx);
                const real c36 = CELL_COEFF_ARTM_BL (cc36, z, x, y, dimmz, dimmx);
                const real c44 = CELL_COEFF_BL      (cc44, z, x, y, dimmz, dimmx);
                const real c45 = CELL_COEFF_ARTM_BL (cc45, z, x, y, dimmz, dimmx);
                const real c46 = CELL_COEFF_ARTM_BL (cc46, z, x, y, dimmz, dimmx);
                const real c55 = CELL_COEFF_BL      (cc55, z, x, y, dimmz, dimmx);
                const real c56 = CELL_COEFF_ARTM_BL (cc56, z, x, y, dimmz, dimmx);
                const real c66 = CELL_COEFF_BL      (cc66, z, x, y, dimmz, dimmx);

                const real u_x = STENCILX (vxu, SX, dxi, z, x, y, dimmz, dimmx);
                const real v_x = STENCILX (vxv, SX, dxi, z, x, y, dimmz, dimmx);
                const real w_x = STENCILX (vxw, SX, dxi, z, x, y, dimmz, dimmx);

                const real u_y = STENCILY (vyu, SY, dyi, z, x, y, dimmz, dimmx);
                const real v_y = STENCILY (vyv, SY, dyi, z, x, y, dimmz, dimmx);
                const real w_y = STENCILY (vyw, SY, dyi, z, x, y, dimmz, dimmx);

                const real u_z = STENCILZ (vzu, SZ, dzi, z, x, y, dimmz, dimmx);
                const real v_z = STENCILZ (vzv, SZ, dzi, z, x, y, dimmz, dimmx);
                const real w_z = STENCILZ (vzw, SZ, dzi, z, x, y, dimmz, dimmx);

                STRESS_UPDATE (sxxptr,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx);
                STRESS_UPDATE (syyptr,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx);
                STRESS_UPDATE (szzptr,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx);
                STRESS_UPDATE (syzptr,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx);
                STRESS_UPDATE (sxzptr,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx);
                STRESS_UPDATE (sxyptr,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx);
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
                       const integer dimmy,
                       const phase_t phase )
{
#if defined(__INTEL_COMPILER)
    #pragma forceinline recursive
#endif
    {
        compute_component_scell_BR ( s, v.tr, v.bl, v.br, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, back_offset, dimmz, dimmx, dimmy, phase);
        compute_component_scell_BL ( s, v.tl, v.br, v.bl, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, forw_offset, dimmz, dimmx, dimmy, phase);
        compute_component_scell_TR ( s, v.br, v.tl, v.tr, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, forw_offset, forw_offset, dimmz, dimmx, dimmy, phase);
        compute_component_scell_TL ( s, v.bl, v.tr, v.tl, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, back_offset, back_offset, dimmz, dimmx, dimmy, phase);
    }
};

