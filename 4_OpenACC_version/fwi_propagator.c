#include "fwi_propagator.h"

inline
integer IDX (const integer z, 
             const integer x, 
             const integer y, 
             const integer dimmz, 
             const integer dimmx)
{
    return (y*dimmx*dimmz) + (x*dimmz) + (z);
};


real stencil_Z ( const integer off,
                 const real* restrict ptr,
                 const real    dzi,
                 const integer z,
                 const integer x,
                 const integer y,
                 const int dimmz,
                 const int dimmx)
{
    return  ((C0 * ( ptr[IDX(z  +off,x,y,dimmz,dimmx)] - ptr[IDX(z-1+off,x,y,dimmz,dimmx)]) +
              C1 * ( ptr[IDX(z+1+off,x,y,dimmz,dimmx)] - ptr[IDX(z-2+off,x,y,dimmz,dimmx)]) +
              C2 * ( ptr[IDX(z+2+off,x,y,dimmz,dimmx)] - ptr[IDX(z-3+off,x,y,dimmz,dimmx)]) +
              C3 * ( ptr[IDX(z+3+off,x,y,dimmz,dimmx)] - ptr[IDX(z-4+off,x,y,dimmz,dimmx)])) * dzi );
};

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
                                 const phase_t        phase)
{
    const integer start   = ((nzf-nz0) + 2*HALO) * ((nxf-nx0) + 2*HALO) * (ny0 - HALO);
    const integer end     = ((nzf-nz0) + 2*HALO) * ((nxf-nx0) + 2*HALO) * (nyf + HALO);
    const integer nelems  = end - start;

    #pragma acc kernels copyin(szptr[start:nelems], sxptr[start:nelems], syptr[start:nelems], rho[start:nelems]) \
                        copyin(vptr[start:nelems]) \
                        async(phase) wait(H2D)
    {
    #pragma acc loop independent 
    for(integer y=ny0; y < nyf; y++)
    {
        #pragma acc loop independent device_type(nvidia) gang worker(4)
        for(integer x=nx0; x < nxf; x++)
        {
            #pragma acc loop independent device_type(nvidia) gang vector(32)
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
    } /* end acc kernels */
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
                                 const phase_t        phase)
{
    const integer start  = ((nzf-nz0) + 2*HALO) * ((nxf-nx0) + 2*HALO) * (ny0 - HALO);
    const integer end    = ((nzf-nz0) + 2*HALO) * ((nxf-nx0) + 2*HALO) * (nyf + HALO);
    const integer nelems = end - start;

    #pragma acc kernels copyin(szptr[start:nelems], sxptr[start:nelems], syptr[start:nelems], rho[start:nelems]) \
                        copyin(vptr[start:nelems]) \
                        async(phase) wait(H2D)
    {
    #pragma acc loop independent 
    for(integer y=ny0; y < nyf; y++)
    {
        #pragma acc loop independent device_type(nvidia) gang worker(4)
        for(integer x=nx0; x < nxf; x++)
        {
            #pragma acc loop independent device_type(nvidia) gang vector(32)
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
    } /* end acc kernels */
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
                                 const phase_t        phase)
{
    const integer start  = ((nzf-nz0) + 2*HALO) * ((nxf-nx0) + 2*HALO) * (ny0 - HALO);
    const integer end    = ((nzf-nz0) + 2*HALO) * ((nxf-nx0) + 2*HALO) * (nyf + HALO);
    const integer nelems = end - start;

    #pragma acc kernels copyin(szptr[start:nelems], sxptr[start:nelems], syptr[start:nelems], rho[start:nelems]) \
                        copyin(vptr[start:nelems]) \
                        async(phase) wait(H2D)
    {
    #pragma acc loop independent
    for(integer y=ny0; y < nyf; y++)
    {
        #pragma acc loop independent device_type(nvidia) gang worker(4)
        for(integer x=nx0; x < nxf; x++)
        {
            #pragma acc loop independent device_type(nvidia) gang vector(32)
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
    } /* end acc kernels */
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
                                 const phase_t        phase)
{
    const integer start  = ((nzf-nz0) + 2*HALO) * ((nxf-nx0) + 2*HALO) * (ny0 - HALO);
    const integer end    = ((nzf-nz0) + 2*HALO) * ((nxf-nx0) + 2*HALO) * (nyf + HALO);
    const integer nelems = end - start;

    #pragma acc kernels copyin(szptr[start:nelems], sxptr[start:nelems], syptr[start:nelems], rho[start:nelems]) \
                        copyin(vptr[start:nelems]) \
                        async(phase) wait(H2D)
    {
    #pragma acc loop independent 
    for(integer y=ny0; y < nyf; y++)
    {
        #pragma acc loop independent device_type(nvidia) gang worker(4)
        for(integer x=nx0; x < nxf; x++)
        {
            #pragma acc loop independent device_type(nvidia) gang vector(32)
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
    } /* end acc loop */
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
                         const phase_t phase)
{
#ifdef DEBUG
    fprintf(stderr, "Integration limits of %s are (z "I"-"I",x "I"-"I",y "I"-"I")\n", __FUNCTION__, nz0,nzf,nx0,nxf,ny0,nyf);
#endif

    compute_component_vcell_TL (v.tl.w, s.bl.zz, s.tr.xz, s.tl.yz, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, back_offset, forw_offset, dimmz, dimmx, phase);
    compute_component_vcell_TR (v.tr.w, s.br.zz, s.tl.xz, s.tr.yz, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, forw_offset, back_offset, dimmz, dimmx, phase);
    compute_component_vcell_BL (v.bl.w, s.tl.zz, s.br.xz, s.bl.yz, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, back_offset, dimmz, dimmx, phase);
    compute_component_vcell_BR (v.br.w, s.tr.zz, s.bl.xz, s.br.yz, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, forw_offset, forw_offset, dimmz, dimmx, phase);
    compute_component_vcell_TL (v.tl.u, s.bl.xz, s.tr.xx, s.tl.xy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, back_offset, forw_offset, dimmz, dimmx, phase);
    compute_component_vcell_TR (v.tr.u, s.br.xz, s.tl.xx, s.tr.xy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, forw_offset, back_offset, dimmz, dimmx, phase);
    compute_component_vcell_BL (v.bl.u, s.tl.xz, s.br.xx, s.bl.xy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, back_offset, dimmz, dimmx, phase);
    compute_component_vcell_BR (v.br.u, s.tr.xz, s.bl.xx, s.br.xy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, forw_offset, forw_offset, dimmz, dimmx, phase);
    compute_component_vcell_TL (v.tl.v, s.bl.yz, s.tr.xy, s.tl.yy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, back_offset, forw_offset, dimmz, dimmx, phase);
    compute_component_vcell_TR (v.tr.v, s.br.yz, s.tl.xy, s.tr.yy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, forw_offset, back_offset, dimmz, dimmx, phase);
    compute_component_vcell_BL (v.bl.v, s.tl.yz, s.br.xy, s.bl.yy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, back_offset, dimmz, dimmx, phase);
    compute_component_vcell_BR (v.br.v, s.tr.yz, s.bl.xy, s.br.yy, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, forw_offset, forw_offset, dimmz, dimmx, phase);
};





/* ------------------------------------------------------------------------------ */
/*                                                                                */
/*                               CALCULO DE TENSIONES                             */
/*                                                                                */
/* ------------------------------------------------------------------------------ */

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
                       const phase_t phase )
{
    compute_component_scell_BR ( s, v.tr, v.bl, v.br, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, back_offset, dimmz, dimmx, phase);
    compute_component_scell_BL ( s, v.tl, v.br, v.bl, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, forw_offset, dimmz, dimmx, phase);
    compute_component_scell_TR ( s, v.br, v.tl, v.tr, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, forw_offset, forw_offset, dimmz, dimmx, phase);
    compute_component_scell_TL ( s, v.bl, v.tr, v.tl, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, back_offset, back_offset, dimmz, dimmx, phase);
};

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

real cell_coeff_TL ( const real* restrict ptr, 
                     const integer z, 
                     const integer x, 
                     const integer y, 
                     const integer dimmz, 
                     const integer dimmx)
{
    return ( 1.0f / (ptr[IDX(z,x,y,dimmz,dimmx)]));
};

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

real cell_coeff_TR ( const real* restrict ptr, 
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

real cell_coeff_ARTM_TL( const real* restrict ptr, 
                         const integer z, 
                         const integer x, 
                         const integer y, 
                         const integer dimmz, 
                         const integer dimmx)
{
    return (1.0f / ptr[IDX(z,x,y,dimmz,dimmx)]);
};

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

    const integer start  = ((nzf-nz0) + 2*HALO) * ((nxf-nx0) + 2*HALO) * (ny0 - HALO);
    const integer end    = ((nzf-nz0) + 2*HALO) * ((nxf-nx0) + 2*HALO) * (nyf + HALO);
    const integer nelems = end - start;

    #pragma acc kernels copyin(sxxptr[start:nelems], syyptr[start:nelems], szzptr[start:nelems], syzptr[start:nelems], sxzptr[start:nelems], sxyptr[start:nelems]) \
                        copyin(vxu[start:nelems], vxv[start:nelems], vxw[start:nelems])  \
                        copyin(vyu[start:nelems], vyv[start:nelems], vyw[start:nelems])  \
                        copyin(vzu[start:nelems], vzv[start:nelems], vzw[start:nelems])  \
                        present(cc11, cc12, cc13, cc14, cc15, cc16)             \
                        present(cc22, cc23, cc24, cc25, cc26)                   \
                        present(cc33, cc34, cc35, cc36)                         \
                        present(cc44, cc45, cc46)                               \
                        present(cc55, cc56)                                     \
                        present(cc66)                                           \
                        async(phase)
    {
    for (integer y = ny0; y < nyf; y++)
    {
        #pragma acc loop device_type(nvidia) gang worker(4)
        for (integer x = nx0; x < nxf; x++)
        {
            #pragma acc loop device_type(nvidia) gang vector(32)
            for (integer z = nz0; z < nzf; z++ )
            {
                const real c11 = cell_coeff_TR      (cc11, z, x, y, dimmz, dimmx);
                const real c12 = cell_coeff_TR      (cc12, z, x, y, dimmz, dimmx);
                const real c13 = cell_coeff_TR      (cc13, z, x, y, dimmz, dimmx);
                const real c14 = cell_coeff_ARTM_TR (cc14, z, x, y, dimmz, dimmx);
                const real c15 = cell_coeff_ARTM_TR (cc15, z, x, y, dimmz, dimmx);
                const real c16 = cell_coeff_ARTM_TR (cc16, z, x, y, dimmz, dimmx);
                const real c22 = cell_coeff_TR      (cc22, z, x, y, dimmz, dimmx);
                const real c23 = cell_coeff_TR      (cc23, z, x, y, dimmz, dimmx);
                const real c24 = cell_coeff_ARTM_TR (cc24, z, x, y, dimmz, dimmx);
                const real c25 = cell_coeff_ARTM_TR (cc25, z, x, y, dimmz, dimmx);
                const real c26 = cell_coeff_ARTM_TR (cc26, z, x, y, dimmz, dimmx);
                const real c33 = cell_coeff_TR      (cc33, z, x, y, dimmz, dimmx);
                const real c34 = cell_coeff_ARTM_TR (cc34, z, x, y, dimmz, dimmx);
                const real c35 = cell_coeff_ARTM_TR (cc35, z, x, y, dimmz, dimmx);
                const real c36 = cell_coeff_ARTM_TR (cc36, z, x, y, dimmz, dimmx);
                const real c44 = cell_coeff_TR      (cc44, z, x, y, dimmz, dimmx);
                const real c45 = cell_coeff_ARTM_TR (cc45, z, x, y, dimmz, dimmx);
                const real c46 = cell_coeff_ARTM_TR (cc46, z, x, y, dimmz, dimmx);
                const real c55 = cell_coeff_TR      (cc55, z, x, y, dimmz, dimmx);
                const real c56 = cell_coeff_ARTM_TR (cc56, z, x, y, dimmz, dimmx);
                const real c66 = cell_coeff_TR      (cc66, z, x, y, dimmz, dimmx);
                
                const real u_x = stencil_X (_SX, vxu, dxi, z, x, y, dimmz, dimmx);
                const real v_x = stencil_X (_SX, vxv, dxi, z, x, y, dimmz, dimmx);
                const real w_x = stencil_X (_SX, vxw, dxi, z, x, y, dimmz, dimmx);
                
                const real u_y = stencil_Y (_SY, vyu, dyi, z, x, y, dimmz, dimmx);
                const real v_y = stencil_Y (_SY, vyv, dyi, z, x, y, dimmz, dimmx);
                const real w_y = stencil_Y (_SY, vyw, dyi, z, x, y, dimmz, dimmx);
                
                const real u_z = stencil_Z (_SZ, vzu, dzi, z, x, y, dimmz, dimmx);
                const real v_z = stencil_Z (_SZ, vzv, dzi, z, x, y, dimmz, dimmx);
                const real w_z = stencil_Z (_SZ, vzw, dzi, z, x, y, dimmz, dimmx);
                
                stress_update (sxxptr,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (syyptr,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (szzptr,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (syzptr,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (sxzptr,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (sxyptr,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            }
        }
    }
    } /* end acc kernels */ 
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

    const integer start  = ((nzf-nz0) + 2*HALO) * ((nxf-nx0) + 2*HALO) * (ny0 - HALO);
    const integer end    = ((nzf-nz0) + 2*HALO) * ((nxf-nx0) + 2*HALO) * (nyf + HALO);
    const integer nelems = end - start;

    #pragma acc kernels copyin(sxxptr[start:nelems], syyptr[start:nelems], szzptr[start:nelems], syzptr[start:nelems], sxzptr[start:nelems], sxyptr[start:nelems]) \
                        copyin(vxu[start:nelems], vxv[start:nelems], vxw[start:nelems])  \
                        copyin(vyu[start:nelems], vyv[start:nelems], vyw[start:nelems])  \
                        copyin(vzu[start:nelems], vzv[start:nelems], vzw[start:nelems])  \
                        present(cc11, cc12, cc13, cc14, cc15, cc16)             \
                        present(cc22, cc23, cc24, cc25, cc26)                   \
                        present(cc33, cc34, cc35, cc36)                         \
                        present(cc44, cc45, cc46)                               \
                        present(cc55, cc56)                                     \
                        present(cc66)                                           \
                        async(phase)
    for (integer y = ny0; y < nyf; y++)
    {
        #pragma acc loop device_type(nvidia) gang worker(4)
        for (integer x = nx0; x < nxf; x++)
        {
            #pragma acc loop device_type(nvidia) gang vector(32)
            for (integer z = nz0; z < nzf; z++ )
            {
                const real c11 = cell_coeff_TL      (cc11, z, x, y, dimmz, dimmx);
                const real c12 = cell_coeff_TL      (cc12, z, x, y, dimmz, dimmx);
                const real c13 = cell_coeff_TL      (cc13, z, x, y, dimmz, dimmx);
                const real c14 = cell_coeff_ARTM_TL (cc14, z, x, y, dimmz, dimmx);
                const real c15 = cell_coeff_ARTM_TL (cc15, z, x, y, dimmz, dimmx);
                const real c16 = cell_coeff_ARTM_TL (cc16, z, x, y, dimmz, dimmx);
                const real c22 = cell_coeff_TL      (cc22, z, x, y, dimmz, dimmx);
                const real c23 = cell_coeff_TL      (cc23, z, x, y, dimmz, dimmx);
                const real c24 = cell_coeff_ARTM_TL (cc24, z, x, y, dimmz, dimmx);
                const real c25 = cell_coeff_ARTM_TL (cc25, z, x, y, dimmz, dimmx);
                const real c26 = cell_coeff_ARTM_TL (cc26, z, x, y, dimmz, dimmx);
                const real c33 = cell_coeff_TL      (cc33, z, x, y, dimmz, dimmx);
                const real c34 = cell_coeff_ARTM_TL (cc34, z, x, y, dimmz, dimmx);
                const real c35 = cell_coeff_ARTM_TL (cc35, z, x, y, dimmz, dimmx);
                const real c36 = cell_coeff_ARTM_TL (cc36, z, x, y, dimmz, dimmx);
                const real c44 = cell_coeff_TL      (cc44, z, x, y, dimmz, dimmx);
                const real c45 = cell_coeff_ARTM_TL (cc45, z, x, y, dimmz, dimmx);
                const real c46 = cell_coeff_ARTM_TL (cc46, z, x, y, dimmz, dimmx);
                const real c55 = cell_coeff_TL      (cc55, z, x, y, dimmz, dimmx);
                const real c56 = cell_coeff_ARTM_TL (cc56, z, x, y, dimmz, dimmx);
                const real c66 = cell_coeff_TL      (cc66, z, x, y, dimmz, dimmx);
                
                const real u_x = stencil_X (_SX, vxu, dxi, z, x, y, dimmz, dimmx);
                const real v_x = stencil_X (_SX, vxv, dxi, z, x, y, dimmz, dimmx);
                const real w_x = stencil_X (_SX, vxw, dxi, z, x, y, dimmz, dimmx);
                
                const real u_y = stencil_Y (_SY, vyu, dyi, z, x, y, dimmz, dimmx);
                const real v_y = stencil_Y (_SY, vyv, dyi, z, x, y, dimmz, dimmx);
                const real w_y = stencil_Y (_SY, vyw, dyi, z, x, y, dimmz, dimmx);
                
                const real u_z = stencil_Z (_SZ, vzu, dzi, z, x, y, dimmz, dimmx);
                const real v_z = stencil_Z (_SZ, vzv, dzi, z, x, y, dimmz, dimmx);
                const real w_z = stencil_Z (_SZ, vzw, dzi, z, x, y, dimmz, dimmx);
                
                stress_update (sxxptr,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (syyptr,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (szzptr,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (syzptr,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (sxzptr,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (sxyptr,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
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

    const integer start  = ((nzf-nz0) + 2*HALO) * ((nxf-nx0) + 2*HALO) * (ny0 - HALO);
    const integer end    = ((nzf-nz0) + 2*HALO) * ((nxf-nx0) + 2*HALO) * (nyf + HALO);
    const integer nelems = end - start;

    #pragma acc kernels copyin(sxxptr[start:nelems], syyptr[start:nelems], szzptr[start:nelems], syzptr[start:nelems], sxzptr[start:nelems], sxyptr[start:nelems]) \
                        copyin(vxu[start:nelems], vxv[start:nelems], vxw[start:nelems])  \
                        copyin(vyu[start:nelems], vyv[start:nelems], vyw[start:nelems])  \
                        copyin(vzu[start:nelems], vzv[start:nelems], vzw[start:nelems])  \
                        present(cc11, cc12, cc13, cc14, cc15, cc16)             \
                        present(cc22, cc23, cc24, cc25, cc26)                   \
                        present(cc33, cc34, cc35, cc36)                         \
                        present(cc44, cc45, cc46)                               \
                        present(cc55, cc56)                                     \
                        present(cc66)                                           \
                        async(phase)
    for (integer y = ny0; y < nyf; y++)
    {
        #pragma acc loop device_type(nvidia) gang worker(4)
        for (integer x = nx0; x < nxf; x++)
        {
            #pragma acc loop device_type(nvidia) gang vector(32)
            for (integer z = nz0; z < nzf; z++ )
            {
                const real c11 = cell_coeff_BR      (cc11, z, x, y, dimmz, dimmx);
                const real c12 = cell_coeff_BR      (cc12, z, x, y, dimmz, dimmx);
                const real c13 = cell_coeff_BR      (cc13, z, x, y, dimmz, dimmx);
                const real c22 = cell_coeff_BR      (cc22, z, x, y, dimmz, dimmx);
                const real c23 = cell_coeff_BR      (cc23, z, x, y, dimmz, dimmx);
                const real c33 = cell_coeff_BR      (cc33, z, x, y, dimmz, dimmx);
                const real c44 = cell_coeff_BR      (cc44, z, x, y, dimmz, dimmx);
                const real c55 = cell_coeff_BR      (cc55, z, x, y, dimmz, dimmx);
                const real c66 = cell_coeff_BR      (cc66, z, x, y, dimmz, dimmx);
                
                const real c14 = cell_coeff_ARTM_BR (cc14, z, x, y, dimmz, dimmx);
                const real c15 = cell_coeff_ARTM_BR (cc15, z, x, y, dimmz, dimmx);
                const real c16 = cell_coeff_ARTM_BR (cc16, z, x, y, dimmz, dimmx);
                const real c24 = cell_coeff_ARTM_BR (cc24, z, x, y, dimmz, dimmx);
                const real c25 = cell_coeff_ARTM_BR (cc25, z, x, y, dimmz, dimmx);
                const real c26 = cell_coeff_ARTM_BR (cc26, z, x, y, dimmz, dimmx);
                const real c34 = cell_coeff_ARTM_BR (cc34, z, x, y, dimmz, dimmx);
                const real c35 = cell_coeff_ARTM_BR (cc35, z, x, y, dimmz, dimmx);
                const real c36 = cell_coeff_ARTM_BR (cc36, z, x, y, dimmz, dimmx);
                const real c45 = cell_coeff_ARTM_BR (cc45, z, x, y, dimmz, dimmx);
                const real c46 = cell_coeff_ARTM_BR (cc46, z, x, y, dimmz, dimmx);
                const real c56 = cell_coeff_ARTM_BR (cc56, z, x, y, dimmz, dimmx);
                
                const real u_x = stencil_X (_SX, vxu, dxi, z, x, y, dimmz, dimmx);
                const real v_x = stencil_X (_SX, vxv, dxi, z, x, y, dimmz, dimmx);
                const real w_x = stencil_X (_SX, vxw, dxi, z, x, y, dimmz, dimmx);
                
                const real u_y = stencil_Y (_SY, vyu, dyi, z, x, y, dimmz, dimmx);
                const real v_y = stencil_Y (_SY, vyv, dyi, z, x, y, dimmz, dimmx);
                const real w_y = stencil_Y (_SY, vyw, dyi, z, x, y, dimmz, dimmx);
                
                const real u_z = stencil_Z (_SZ, vzu, dzi, z, x, y, dimmz, dimmx);
                const real v_z = stencil_Z (_SZ, vzv, dzi, z, x, y, dimmz, dimmx);
                const real w_z = stencil_Z (_SZ, vzw, dzi, z, x, y, dimmz, dimmx);
                
                stress_update (sxxptr,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (syyptr,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (szzptr,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (syzptr,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (sxzptr,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (sxyptr,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
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

    const integer start  = ((nzf-nz0) + 2*HALO) * ((nxf-nx0) + 2*HALO) * (ny0 - HALO);
    const integer end    = ((nzf-nz0) + 2*HALO) * ((nxf-nx0) + 2*HALO) * (nyf + HALO);
    const integer nelems = end - start;

    #pragma acc kernels copyin(sxxptr[start:nelems], syyptr[start:nelems], szzptr[start:nelems], syzptr[start:nelems], sxzptr[start:nelems], sxyptr[start:nelems]) \
                        copyin(vxu[start:nelems], vxv[start:nelems], vxw[start:nelems])  \
                        copyin(vyu[start:nelems], vyv[start:nelems], vyw[start:nelems])  \
                        copyin(vzu[start:nelems], vzv[start:nelems], vzw[start:nelems])  \
                        present(cc11, cc12, cc13, cc14, cc15, cc16)             \
                        present(cc22, cc23, cc24, cc25, cc26)                   \
                        present(cc33, cc34, cc35, cc36)                         \
                        present(cc44, cc45, cc46)                               \
                        present(cc55, cc56)                                     \
                        present(cc66)                                           \
                        async(phase)
    for (integer y = ny0; y < nyf; y++)
    {
        #pragma acc loop device_type(nvidia) gang worker(4)
        for (integer x = nx0; x < nxf; x++)
        {
            #pragma acc loop device_type(nvidia) gang vector(32)
            for (integer z = nz0; z < nzf; z++ )
            {
                const real c11 = cell_coeff_BL      (cc11, z, x, y, dimmz, dimmx);
                const real c12 = cell_coeff_BL      (cc12, z, x, y, dimmz, dimmx);
                const real c13 = cell_coeff_BL      (cc13, z, x, y, dimmz, dimmx);
                const real c14 = cell_coeff_ARTM_BL (cc14, z, x, y, dimmz, dimmx);
                const real c15 = cell_coeff_ARTM_BL (cc15, z, x, y, dimmz, dimmx);
                const real c16 = cell_coeff_ARTM_BL (cc16, z, x, y, dimmz, dimmx);
                const real c22 = cell_coeff_BL      (cc22, z, x, y, dimmz, dimmx);
                const real c23 = cell_coeff_BL      (cc23, z, x, y, dimmz, dimmx);
                const real c24 = cell_coeff_ARTM_BL (cc24, z, x, y, dimmz, dimmx);
                const real c25 = cell_coeff_ARTM_BL (cc25, z, x, y, dimmz, dimmx);
                const real c26 = cell_coeff_ARTM_BL (cc26, z, x, y, dimmz, dimmx);
                const real c33 = cell_coeff_BL      (cc33, z, x, y, dimmz, dimmx);
                const real c34 = cell_coeff_ARTM_BL (cc34, z, x, y, dimmz, dimmx);
                const real c35 = cell_coeff_ARTM_BL (cc35, z, x, y, dimmz, dimmx);
                const real c36 = cell_coeff_ARTM_BL (cc36, z, x, y, dimmz, dimmx);
                const real c44 = cell_coeff_BL      (cc44, z, x, y, dimmz, dimmx);
                const real c45 = cell_coeff_ARTM_BL (cc45, z, x, y, dimmz, dimmx);
                const real c46 = cell_coeff_ARTM_BL (cc46, z, x, y, dimmz, dimmx);
                const real c55 = cell_coeff_BL      (cc55, z, x, y, dimmz, dimmx);
                const real c56 = cell_coeff_ARTM_BL (cc56, z, x, y, dimmz, dimmx);
                const real c66 = cell_coeff_BL      (cc66, z, x, y, dimmz, dimmx);
                
                const real u_x = stencil_X (_SX, vxu, dxi, z, x, y, dimmz, dimmx);
                const real v_x = stencil_X (_SX, vxv, dxi, z, x, y, dimmz, dimmx);
                const real w_x = stencil_X (_SX, vxw, dxi, z, x, y, dimmz, dimmx);
                
                const real u_y = stencil_Y (_SY, vyu, dyi, z, x, y, dimmz, dimmx);
                const real v_y = stencil_Y (_SY, vyv, dyi, z, x, y, dimmz, dimmx);
                const real w_y = stencil_Y (_SY, vyw, dyi, z, x, y, dimmz, dimmx);
                
                const real u_z = stencil_Z (_SZ, vzu, dzi, z, x, y, dimmz, dimmx);
                const real v_z = stencil_Z (_SZ, vzv, dzi, z, x, y, dimmz, dimmx);
                const real w_z = stencil_Z (_SZ, vzw, dzi, z, x, y, dimmz, dimmx);
                
                stress_update (sxxptr,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx);
                stress_update (syyptr,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx);
                stress_update (szzptr,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx);
                stress_update (syzptr,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx);
                stress_update (sxzptr,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx);
                stress_update (sxyptr,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx);
            }
        }
    }
};
