#include "fwi_propagator.h"

integer IDX(const integer z, 
						const integer x, 
						const integer y, 
						const integer dimmz, 
						const integer dimmx)
{
	return (y*dimmx*dimmz) + (x*dimmz) + (z);
};

real stencil_Z (  const offset_t off,
                  real* restrict ptr,
                  const real    dzi,
                  const integer z,
                  const integer x,
                  const integer y,
                  const integer dimmz,
                  const integer dimmx)
{
  return  ((C0 * ( 	ptr[IDX(z  +off,x,y,dimmz,dimmx)]  - 
										ptr[IDX(z-1+off,x,y,dimmz,dimmx)]) +
            C1 * ( 	ptr[IDX(z+1+off,x,y,dimmz,dimmx)]  - 
										ptr[IDX(z-2+off,x,y,dimmz,dimmx)]) +
            C2 * ( 	ptr[IDX(z+2+off,x,y,dimmz,dimmx)]  - 
										ptr[IDX(z-3+off,x,y,dimmz,dimmx)]) +
            C3 * ( 	ptr[IDX(z+3+off,x,y,dimmz,dimmx)]  - 
										ptr[IDX(z-4+off,x,y,dimmz,dimmx)])) * dzi );
};

real stencil_X( const offset_t off,
                real* restrict ptr,
                const real dxi,
                const integer z,
                const integer x,
                const integer y,
                const integer dimmz,
                const integer dimmx)
{
	return ( (C0 * ( 	ptr[IDX(z,x  +off,y,dimmz,dimmx)]  - 
										ptr[IDX(z,x-1+off,y,dimmz,dimmx)]) +
            C1 * ( 	ptr[IDX(z,x+1+off,y,dimmz,dimmx)]  - 
										ptr[IDX(z,x-2+off,y,dimmz,dimmx)]) +
            C2 * ( 	ptr[IDX(z,x+2+off,y,dimmz,dimmx)]  - 
										ptr[IDX(z,x-3+off,y,dimmz,dimmx)]) +
            C3 * ( 	ptr[IDX(z,x+3+off,y,dimmz,dimmx)]  - 
										ptr[IDX(z,x-4+off,y,dimmz,dimmx)])) * dxi );
};

real stencil_Y( const offset_t off,
                real* restrict ptr,
                const real dyi,
                const integer z,
                const integer x,
                const integer y,
                const integer dimmz,
                const integer dimmx)
{
  return ( (C0 * ( ptr[IDX(z,x,y  +off,dimmz,dimmx)]  -
									 ptr[IDX(z,x,y-1+off,dimmz,dimmx)]) +
            C1 * ( ptr[IDX(z,x,y+1+off,dimmz,dimmx)]  -
						       ptr[IDX(z,x,y-2+off,dimmz,dimmx)]) +
            C2 * ( ptr[IDX(z,x,y+2+off,dimmz,dimmx)]  -
						       ptr[IDX(z,x,y-3+off,dimmz,dimmx)]) +
            C3 * ( ptr[IDX(z,x,y+3+off,dimmz,dimmx)]  - 
						       ptr[IDX(z,x,y-4+off,dimmz,dimmx)])) * dyi );
};

/* -------------------------------------------------------------------- */
/*                     KERNELS FOR VELOCITY                             */
/* -------------------------------------------------------------------- */


real rho_BL ( real* restrict rho,
              const integer z,
              const integer x,
              const integer y,
              const integer dimmz,
              const integer dimmx)
{
    return (2.0f / (rho[IDX(z,x,y,dimmz,dimmx)] + 
										rho[IDX(z+1,x,y,dimmz,dimmx)]));
};

real rho_TR ( real* restrict rho,
              const integer z,
              const integer x,
              const integer y,
              const integer dimmz,
              const integer dimmx)
{
    return (2.0f/ (rho[IDX(z,x,y,dimmz,dimmx)] +  
									 rho[IDX(z,x+1,y,dimmz,dimmx)]) );
};

real rho_BR ( real* restrict rho,
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

real rho_TL ( real* restrict rho,
              const integer z,
              const integer x,
              const integer y,
              const integer dimmz,
              const integer dimmx)
{
    return (2.0f/ (rho[IDX(z,x,y,dimmz,dimmx)] + 
					         rho[IDX(z,x,y+1,dimmz,dimmx)]));
};

void compute_component_vcell_TL (real* restrict vptr,
                                 real* restrict szptr,
                                 real* restrict sxptr,
                                 real* restrict syptr,
                                 real* restrict rho,
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
	__assume( nz0 % HALO == 0);
	__assume( nzf % HALO == 0);

	real* restrict _vptr  __attribute__ ((aligned (64))) = vptr ;
	real* restrict _szptr __attribute__ ((aligned (64))) = szptr;
	real* restrict _sxptr __attribute__ ((aligned (64))) = sxptr;
	real* restrict _syptr __attribute__ ((aligned (64))) = syptr;

	#pragma omp parallel for
  for(integer y=ny0; y < nyf; y++)
    for(integer x=nx0; x < nxf; x++)
		{
			#pragma simd
      for(integer z=nz0; z < nzf; z++)
      {
        const real lrho = rho_TL(rho, z, x, y, dimmz, dimmx);

        const real stx  = stencil_X( _SX, _sxptr, dxi, z, x, y, dimmz, dimmx);
        const real sty  = stencil_Y( _SY, _syptr, dyi, z, x, y, dimmz, dimmx);
        const real stz  = stencil_Z( _SZ, _szptr, dzi, z, x, y, dimmz, dimmx);

        _vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
      }
		}
};

void compute_component_vcell_TR (real* restrict vptr,
                                 real* restrict szptr,
                                 real* restrict sxptr,
                                 real* restrict syptr,
                                 real* restrict rho,
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
	__assume( nz0 % HALO == 0);
	__assume( nzf % HALO == 0);

	real* restrict _vptr  __attribute__ ((aligned (64))) = vptr ;
	real* restrict _szptr __attribute__ ((aligned (64))) = szptr;
	real* restrict _sxptr __attribute__ ((aligned (64))) = sxptr;
	real* restrict _syptr __attribute__ ((aligned (64))) = syptr;

	#pragma omp parallel for
  for(integer y=ny0; y < nyf; y++)
    for(integer x=nx0; x < nxf; x++)
		{
			#pragma simd
			for(integer z=nz0; z < nzf; z++)
        {
          const real lrho = rho_TR(rho, z, x, y, dimmz, dimmx);

          const real stx  = stencil_X( _SX, _sxptr, dxi, z, x, y, dimmz, dimmx);
          const real sty  = stencil_Y( _SY, _syptr, dyi, z, x, y, dimmz, dimmx);
          const real stz  = stencil_Z( _SZ, _szptr, dzi, z, x, y, dimmz, dimmx);

          _vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
        }
			}
};

void compute_component_vcell_BR (real* restrict  vptr,
                                 real* restrict  szptr,
                                 real* restrict  sxptr,
                                 real* restrict  syptr,
                                 real* restrict  rho,
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
	__assume( nz0 % HALO == 0);
	__assume( nzf % HALO == 0);
	
	real* restrict _vptr  __attribute__ ((aligned (64))) = vptr ;
	real* restrict _szptr __attribute__ ((aligned (64))) = szptr;
	real* restrict _sxptr __attribute__ ((aligned (64))) = sxptr;
	real* restrict _syptr __attribute__ ((aligned (64))) = syptr;


	#pragma omp parallel for
  for(integer y=ny0; y < nyf; y++)
  	for(integer x=nx0; x < nxf; x++)
		{
			#pragma simd
      for(integer z=nz0; z < nzf; z++)
      {
          const real lrho = rho_BR(rho, z, x, y, dimmz, dimmx);

          const real stx  = stencil_X( _SX, _sxptr, dxi, z, x, y, dimmz, dimmx );
          const real sty  = stencil_Y( _SY, _syptr, dyi, z, x, y, dimmz, dimmx );
          const real stz  = stencil_Z( _SZ, _szptr, dzi, z, x, y, dimmz, dimmx );

          _vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
      }
		}
};

void compute_component_vcell_BL (real* restrict  vptr,
                                 real* restrict  szptr,
                                 real* restrict  sxptr,
                                 real* restrict  syptr,
                                 real* restrict  rho,
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
	__assume( nz0 % HALO == 0);
	__assume( nzf % HALO == 0);
	
	real* restrict _vptr  __attribute__ ((aligned (64))) = vptr ;
	real* restrict _szptr __attribute__ ((aligned (64))) = szptr;
	real* restrict _sxptr __attribute__ ((aligned (64))) = sxptr;
	real* restrict _syptr __attribute__ ((aligned (64))) = syptr;


	#pragma omp parallel for
  for(integer y=ny0; y < nyf; y++)
  	for(integer x=nx0; x < nxf; x++)
		{
			#pragma simd
      for(integer z=nz0; z < nzf; z++)
      {
        const real lrho = rho_TL(rho, z, x, y, dimmz, dimmx);

        const real stx  = stencil_X( _SX, _sxptr, dxi, z, x, y, dimmz, dimmx);
        const real sty  = stencil_Y( _SY, _syptr, dyi, z, x, y, dimmz, dimmx);
        const real stz  = stencil_Z( _SZ, _szptr, dzi, z, x, y, dimmz, dimmx);

        _vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
      }
		}
};

void velocity_propagator(v_t       v,
                         s_t       s,
                         coeff_t   coeffs,
                         real      *rho,
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
                         const integer   dimmz,
                         const integer   dimmx)
{
#ifdef DEBUG
	fprintf(stderr, "Integration limits of %s are (z "I"-"I",x "I"-"I",y "I"-"I")\n", __FUNCTION__, nz0,nzf,nx0,nxf,ny0,nyf);
#endif

	#pragma forceinline recursive
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
};





/* ------------------------------------------------------------------------------ */
/*                                                                                */
/*                               CALCULO DE TENSIONES                             */
/*                                                                                */
/* ------------------------------------------------------------------------------ */

void stress_update( real* restrict sptr,
                   const real          c1,
                   const real          c2,
                   const real          c3,
                   const real          c4,
                   const real          c5,
                   const real          c6,
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
  sptr[IDX(z,x,y,dimmz,dimmx)] += dt * c1 * u_x;
  sptr[IDX(z,x,y,dimmz,dimmx)] += dt * c2 * v_y;
  sptr[IDX(z,x,y,dimmz,dimmx)] += dt * c3 * w_z;
  sptr[IDX(z,x,y,dimmz,dimmx)] += dt * c4 * (w_y + v_z);
  sptr[IDX(z,x,y,dimmz,dimmx)] += dt * c5 * (w_x + u_z);
  sptr[IDX(z,x,y,dimmz,dimmx)] += dt * c6 * (v_x + u_y);
};

void stress_propagator( s_t          s,
                       v_t           v,
                       coeff_t       coeffs,
                       real          *rho,
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
                       const integer dimmx )
{
	#pragma forceinline recursive
  {
    compute_component_scell_BR ( s, v.tr, v.bl, v.br, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, back_offset, dimmz, dimmx);
    compute_component_scell_BL ( s, v.tl, v.br, v.bl, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, forw_offset, back_offset, forw_offset, dimmz, dimmx);
    compute_component_scell_TR ( s, v.br, v.tl, v.tr, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, forw_offset, forw_offset, dimmz, dimmx);
    compute_component_scell_TL ( s, v.bl, v.tr, v.tl, coeffs, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, back_offset, back_offset, back_offset, dimmz, dimmx);
  }
};

real cell_coeff_BR ( real* restrict ptr, const integer z, const integer x, const integer y, const integer dimmz, const integer dimmx)
{
    return ( 1.0f / ( 2.5f  *(ptr[IDX(z  , x  ,y,dimmz,dimmx)] +
                            ptr[IDX(z  , x+1,y,dimmz,dimmx)] +
                            ptr[IDX(z+1, x  ,y,dimmz,dimmx)] +
                            ptr[IDX(z+1, x+1,y,dimmz,dimmx)])) );
};

real cell_coeff_TL ( real* restrict ptr, const integer z, const integer x, const integer y, const integer dimmz, const integer dimmx)
{
    return ( 1.0f / (ptr[IDX(z,x,y,dimmz,dimmx)]));
};

real cell_coeff_BL ( real* restrict ptr, const integer z, const integer x, const integer y, const integer dimmz, const integer dimmx)
{
	return ( 1.0f / ( 2.5f *(ptr[IDX(z  ,x,y  ,dimmz,dimmx)] +
                         ptr[IDX(z  ,x,y+1,dimmz,dimmx)] +
                         ptr[IDX(z+1,x,y  ,dimmz,dimmx)] +
                         ptr[IDX(z+1,x,y+1,dimmz,dimmx)])) );
};

real cell_coeff_TR ( real* restrict ptr, const integer z, const integer x, const integer y, const integer dimmz, const integer dimmx)
{
    return ( 1.0f / ( 2.5f *(ptr[IDX(z  , x  , y  ,dimmz,dimmx)] +
                            ptr[IDX(z  , x+1, y  ,dimmz,dimmx)] +
                            ptr[IDX(z  , x  , y+1,dimmz,dimmx)] +
                            ptr[IDX(z  , x+1, y+1,dimmz,dimmx)])));
};

real cell_coeff_ARTM_BR( real* restrict ptr, const integer z, const integer x, const integer y, const integer dimmz, const integer dimmx)
{
    return ((1.0f / ptr[IDX(z  ,x  ,y,dimmz,dimmx )]  +
             1.0f / ptr[IDX(z  ,x+1,y,dimmz,dimmx )]  +
             1.0f / ptr[IDX(z+1,x  ,y,dimmz,dimmx )]  +
             1.0f / ptr[IDX(z+1,x+1,y,dimmz,dimmx )]) * 0.25f);
};

real cell_coeff_ARTM_TL( real* restrict ptr, const integer z, const integer x, const integer y, const integer dimmz, const integer dimmx)
{
    return (1.0f / ptr[IDX(z,x,y,dimmz,dimmx)]);
};

real cell_coeff_ARTM_BL( real* restrict ptr, const integer z, const integer x, const integer y, const integer dimmz, const integer dimmx)
{
    return ((1.0f / ptr[IDX(z  ,x,y  ,dimmz,dimmx)]  +
             1.0f / ptr[IDX(z  ,x,y+1,dimmz,dimmx)]  +
             1.0f / ptr[IDX(z+1,x,y  ,dimmz,dimmx)]  +
             1.0f / ptr[IDX(z+1,x,y+1,dimmz,dimmx)]) * 0.25f);
};

real cell_coeff_ARTM_TR( real* restrict ptr, const integer z, const integer x, const integer y, const integer dimmz, const integer dimmx)
{
    return ((1.0f / ptr[IDX(z,x  ,y  ,dimmz,dimmx)]  +
             1.0f / ptr[IDX(z,x+1,y  ,dimmz,dimmx)]  +
             1.0f / ptr[IDX(z,x  ,y+1,dimmz,dimmx)]  +
             1.0f / ptr[IDX(z,x+1,y+1,dimmz,dimmx)]) * 0.25f);
};

void compute_component_scell_TR ( s_t             s,
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
	__assume( nz0 % HALO == 0);
	__assume( nzf % HALO == 0);

  real* restrict sxxptr __attribute__ ((aligned (64))) = s.tr.xx;
  real* restrict syyptr __attribute__ ((aligned (64))) = s.tr.yy;
  real* restrict szzptr __attribute__ ((aligned (64))) = s.tr.zz;
  real* restrict syzptr __attribute__ ((aligned (64))) = s.tr.yz;
  real* restrict sxzptr __attribute__ ((aligned (64))) = s.tr.xz;
  real* restrict sxyptr __attribute__ ((aligned (64))) = s.tr.xy;

  real* restrict vxu    __attribute__ ((aligned (64))) = vnode_x.u;
  real* restrict vxv    __attribute__ ((aligned (64))) = vnode_x.v;
  real* restrict vxw    __attribute__ ((aligned (64))) = vnode_x.w;

  real* restrict vyu    __attribute__ ((aligned (64))) = vnode_y.u;
  real* restrict vyv    __attribute__ ((aligned (64))) = vnode_y.v;
  real* restrict vyw    __attribute__ ((aligned (64))) = vnode_y.w;
  real* restrict vzu    __attribute__ ((aligned (64))) = vnode_z.u;
  real* restrict vzv    __attribute__ ((aligned (64))) = vnode_z.v;
  real* restrict vzw    __attribute__ ((aligned (64))) = vnode_z.w;

	#pragma omp parallel for
  for (integer y = ny0; y < nyf; y++)
  	for (integer x = nx0; x < nxf; x++)
		{
			#pragma simd
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
};

void compute_component_scell_TL ( s_t             s,
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
	__assume( nz0 % HALO == 0);
	__assume( nzf % HALO == 0);

  real* restrict sxxptr __attribute__ ((aligned (64))) = s.tl.xx;
  real* restrict syyptr __attribute__ ((aligned (64))) = s.tl.yy;
  real* restrict szzptr __attribute__ ((aligned (64))) = s.tl.zz;
  real* restrict syzptr __attribute__ ((aligned (64))) = s.tl.yz;
  real* restrict sxzptr __attribute__ ((aligned (64))) = s.tl.xz;
  real* restrict sxyptr __attribute__ ((aligned (64))) = s.tl.xy;

  real* restrict vxu    __attribute__ ((aligned (64))) = vnode_x.u;
  real* restrict vxv    __attribute__ ((aligned (64))) = vnode_x.v;
  real* restrict vxw    __attribute__ ((aligned (64))) = vnode_x.w;

  real* restrict vyu    __attribute__ ((aligned (64))) = vnode_y.u;
  real* restrict vyv    __attribute__ ((aligned (64))) = vnode_y.v;
  real* restrict vyw    __attribute__ ((aligned (64))) = vnode_y.w;

  real* restrict vzu    __attribute__ ((aligned (64))) = vnode_z.u;
  real* restrict vzv    __attribute__ ((aligned (64))) = vnode_z.v;
  real* restrict vzw    __attribute__ ((aligned (64))) = vnode_z.w;

	#pragma omp parallel for
  for (integer y = ny0; y < nyf; y++)
  	for (integer x = nx0; x < nxf; x++)
		{
			#pragma simd
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
};


void compute_component_scell_BR ( s_t             s,
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
	__assume( nz0 % HALO == 0);
	__assume( nzf % HALO == 0);

  real* restrict sxxptr __attribute__ ((aligned (64))) = s.br.xx;
  real* restrict syyptr __attribute__ ((aligned (64))) = s.br.yy;
  real* restrict szzptr __attribute__ ((aligned (64))) = s.br.zz;
  real* restrict syzptr __attribute__ ((aligned (64))) = s.br.yz;
  real* restrict sxzptr __attribute__ ((aligned (64))) = s.br.xz;
  real* restrict sxyptr __attribute__ ((aligned (64))) = s.br.xy;

  real* restrict vxu    __attribute__ ((aligned (64))) = vnode_x.u;
  real* restrict vxv    __attribute__ ((aligned (64))) = vnode_x.v;
  real* restrict vxw    __attribute__ ((aligned (64))) = vnode_x.w;

  real* restrict vyu    __attribute__ ((aligned (64))) = vnode_y.u;
  real* restrict vyv    __attribute__ ((aligned (64))) = vnode_y.v;
  real* restrict vyw    __attribute__ ((aligned (64))) = vnode_y.w;

  real* restrict vzu    __attribute__ ((aligned (64))) = vnode_z.u;
  real* restrict vzv    __attribute__ ((aligned (64))) = vnode_z.v;
  real* restrict vzw    __attribute__ ((aligned (64))) = vnode_z.w;

	#pragma omp parallel for
  for (integer y = ny0; y < nyf; y++)
  	for (integer x = nx0; x < nxf; x++)
		{
			#pragma simd
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
};

void compute_component_scell_BL ( s_t             s,
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
	__assume( nz0 % HALO == 0);
	__assume( nzf % HALO == 0);

  real* restrict sxxptr __attribute__ ((aligned (64))) = s.br.xx;
  real* restrict syyptr __attribute__ ((aligned (64))) = s.br.yy;
  real* restrict szzptr __attribute__ ((aligned (64))) = s.br.zz;
  real* restrict syzptr __attribute__ ((aligned (64))) = s.br.yz;
  real* restrict sxzptr __attribute__ ((aligned (64))) = s.br.xz;
  real* restrict sxyptr __attribute__ ((aligned (64))) = s.br.xy;

  real* restrict vxu    __attribute__ ((aligned (64))) = vnode_x.u;
  real* restrict vxv    __attribute__ ((aligned (64))) = vnode_x.v;
  real* restrict vxw    __attribute__ ((aligned (64))) = vnode_x.w;

  real* restrict vyu    __attribute__ ((aligned (64))) = vnode_y.u;
  real* restrict vyv    __attribute__ ((aligned (64))) = vnode_y.v;
  real* restrict vyw    __attribute__ ((aligned (64))) = vnode_y.w;

  real* restrict vzu    __attribute__ ((aligned (64))) = vnode_z.u;
  real* restrict vzv    __attribute__ ((aligned (64))) = vnode_z.v;
  real* restrict vzw    __attribute__ ((aligned (64))) = vnode_z.w;

	#pragma omp parallel for
  for (integer y = ny0; y < nyf; y++)
  	for (integer x = nx0; x < nxf; x++)
		{
			#pragma simd
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
};
