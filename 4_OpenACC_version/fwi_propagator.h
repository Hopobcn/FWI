#ifndef _FWI_PROPAGATOR_H_
#define _FWI_PROPAGATOR_H_

#include "fwi_common.h"

/* stress point structure */
typedef struct {
    real *zz, *xz, *yz, *xx, *xy, *yy;
} point_s_t;

/* velocity point structure */
typedef struct {
    real *u, *v, *w;
} point_v_t;

/* velocity points on a cell */
typedef struct {
    point_v_t tl, tr, bl, br;
} v_t;

/* stress points on a cell */
typedef struct {
    point_s_t tl, tr, bl, br;
} s_t;

/* coefficients for materials */
typedef struct {
    real *c11, *c12, *c13, *c14, *c15, *c16;
    real *c22, *c23, *c24, *c25, *c26;
    real *c33, *c34, *c35, *c36;
    real *c44, *c45, *c46;
    real *c55, *c56;
    real *c66;
} coeff_t;

#define C0 1.2f
#define C1 1.4f
#define C2 1.6f
#define C3 1.8f

#define ASSUMED_DISTANCE 16

typedef enum {back_offset, forw_offset} offset_t;
typedef enum {ONE_A, ONE_B, TWO} phase_t;

#pragma acc routine seq
integer IDX (const integer z, 
             const integer x, 
             const integer y, 
             const integer dimmz, 
             const integer dimmx);

#pragma acc routine seq
real stencil_Z(const integer off,
               const real* restrict ptr,
               const real    dzi,
               const integer z,
               const integer x,
               const integer y,
               const integer dimmz,
               const integer dimmx);

#pragma acc routine seq
real stencil_X(const integer off,
               const real* restrict ptr,
               const real dxi,
               const integer z,
               const integer x,
               const integer y,
               const integer dimmz,
               const integer dimmx);

#pragma acc routine seq
real stencil_Y(const integer off,
               const real* restrict ptr,
               const real dyi,
               const integer z,
               const integer x,
               const integer y,
               const integer dimmz,
               const integer dimmx);


/* ------------------------------------------------------------------------------ */
/*                                                                                */
/*                               CALCULO DE VELOCIDADES                           */
/*                                                                                */
/* ------------------------------------------------------------------------------ */


real rho_BL ( const real* restrict rho,
              const integer z,
              const integer x,
              const integer y,
              const integer dimmz,
              const integer dimmx);

real rho_TR ( const real* restrict rho,
              const integer z,
              const integer x,
              const integer y,
              const integer dimmz,
              const integer dimmx);

real rho_BR ( const real* restrict rho,
              const integer z,
              const integer x,
              const integer y,
              const integer dimmz,
              const integer dimmx);

real rho_TL ( const real* restrict rho,
              const integer z,
              const integer x,
              const integer y,
              const integer dimmz,
              const integer dimmx);

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
                                 const phase_t        phase);


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
                                 const phase_t        phase);


void compute_component_vcell_BR (      real* restrict vptr,
                                 const real* restrict szptr,
                                 const real* restrict sxptr,
                                 const real* restrict syptr,
                                 const real* restrict rho,
                                 const real           dt,
                                 const real           dzi,
                                 const real           dxi,
                                 const real           dyi,
                                 const integer        ny0,
                                 const integer        nyf,
                                 const integer        nx0,
                                 const integer        nxf,
                                 const integer        nz0,
                                 const integer        nzf,
                                 const offset_t       _SZ,
                                 const offset_t       _SX,
                                 const offset_t       _SY,
                                 const integer        dimmz,
                                 const integer        dimmx,
                                 const phase_t        phase);

void compute_component_vcell_BL (      real* restrict vptr,
                                 const real* restrict szptr,
                                 const real* restrict sxptr,
                                 const real* restrict syptr,
                                 const real* restrict rho,
                                 const real           dt,
                                 const real           dzi,
                                 const real           dxi,
                                 const real           dyi,
                                 const integer        ny0,
                                 const integer        nyf,
                                 const integer        nx0,
                                 const integer        nxf,
                                 const integer        nz0,
                                 const integer        nzf,
                                 const offset_t       _SZ,
                                 const offset_t       _SX,
                                 const offset_t       _SY,
                                 const integer        dimmz,
                                 const integer        dimmx,
                                 const phase_t        phase);

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
                         const phase_t phase);





/* ------------------------------------------------------------------------------ */
/*                                                                                */
/*                               CALCULO DE TENSIONES                             */
/*                                                                                */
/* ------------------------------------------------------------------------------ */

#pragma acc routine seq
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
                   const integer  dimmx);

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
                       const phase_t phase );

#pragma acc routine seq
real cell_coeff_BR ( const real* restrict ptr, 
                     const integer z, 
                     const integer x, 
                     const integer y, 
                     const integer dimmz, 
                     const integer dimmx );
#pragma acc routine seq
real cell_coeff_TL ( const real* restrict ptr, 
                     const integer z, 
                     const integer x, 
                     const integer y, 
                     const integer dimmz, 
                     const integer dimmx );
#pragma acc routine seq
real cell_coeff_BL ( const real* restrict ptr, 
                     const integer z, 
                     const integer x, 
                     const integer y, 
                     const integer dimmz, 
                     const integer dimmx );
#pragma acc routine seq
real cell_coeff_TR ( const real* restrict ptr, 
                     const integer z, 
                     const integer x, 
                     const integer y, 
                     const integer dimmz, 
                     const integer dimmx );

#pragma acc routine seq
real cell_coeff_ARTM_BR ( const real* restrict ptr, 
                          const integer z, 
                          const integer x, 
                          const integer y, 
                          const integer dimmz, 
                          const integer dimmx);
#pragma acc routine seq
real cell_coeff_ARTM_TL ( const real* restrict ptr, 
                          const integer z, 
                          const integer x, 
                          const integer y, 
                          const integer dimmz, 
                          const integer dimmx);
#pragma acc routine seq
real cell_coeff_ARTM_BL ( const real* restrict ptr, 
                          const integer z, 
                          const integer x, 
                          const integer y, 
                          const integer dimmz, 
                          const integer dimmx);
#pragma acc routine seq
real cell_coeff_ARTM_TR ( const real* restrict ptr, 
                          const integer z, 
                          const integer x, 
                          const integer y, 
                          const integer dimmz, 
                          const integer dimmx);

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
                                 const integer  dimmx,
                                 const phase_t phase);

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
                                  const integer   dimmz,
                                  const integer   dimmx,
                                  const phase_t   phase);

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
                                  const integer   dimmz,
                                  const integer   dimmx,
                                  const phase_t   phase);

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
                                  const integer   dimmz,
                                  const integer   dimmx,
                                  const phase_t   phase);

#endif /* end of _FWI_PROPAGATOR_H_ definition */
