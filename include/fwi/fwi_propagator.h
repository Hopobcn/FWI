/*
 * =============================================================================
 * Copyright (c) 2016-2018, Barcelona Supercomputing Center (BSC)
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

#ifndef _FWI_PROPAGATOR_H_
#define _FWI_PROPAGATOR_H_

#ifdef __cplusplus
extern "C"
{
#endif

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
typedef enum {ONE_R, ONE_L, TWO, H2D, D2H} phase_t;

#if defined(_OPENACC) 
#pragma acc routine seq
#endif
integer IDX (const integer z, 
             const integer x, 
             const integer y, 
             const integer dimmz, 
             const integer dimmx);

#if defined(_OPENACC)
#pragma acc routine seq
#endif
real stencil_Z(const integer off,
               const real* restrict ptr,
               const real    dzi,
               const integer z,
               const integer x,
               const integer y,
               const integer dimmz,
               const integer dimmx);

#if defined(_OPENACC)
#pragma acc routine seq
#endif
real stencil_X(const integer off,
               const real* restrict ptr,
               const real dxi,
               const integer z,
               const integer x,
               const integer y,
               const integer dimmz,
               const integer dimmx);

#if defined(_OPENACC)
#pragma acc routine seq
#endif
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

#if defined(_OPENACC)
#pragma acc routine seq
#endif
real rho_BL ( const real* restrict rho,
              const integer z,
              const integer x,
              const integer y,
              const integer dimmz,
              const integer dimmx);

#if defined(_OPENACC)
#pragma acc routine seq
#endif
real rho_TR ( const real* restrict rho,
              const integer z,
              const integer x,
              const integer y,
              const integer dimmz,
              const integer dimmx);

#if defined(_OPENACC)
#pragma acc routine seq
#endif
real rho_BR ( const real* restrict rho,
              const integer z,
              const integer x,
              const integer y,
              const integer dimmz,
              const integer dimmx);

#if defined(_OPENACC)
#pragma acc routine seq
#endif
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

/* CUDA STUB */
extern
void compute_component_vcell_TL_cuda ( float* vptr,
                                 const float* szptr,
                                 const float* sxptr,
                                 const float* syptr,
                                 const float* rho,
                                 const float  dt,
                                 const float  dzi,
                                 const float  dxi,
                                 const float  dyi,
                                 const int    nz0,
                                 const int    nzf,
                                 const int    nx0,
                                 const int    nxf,
                                 const int    ny0,
                                 const int    nyf,
                                 const int    SZ,
                                 const int    SX,
                                 const int    SY,
                                 const int    dimmz,
                                 const int    dimmx,
                                 void*        stream);


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
/* CUDA STUB */
extern
void compute_component_vcell_TR_cuda ( float* vptr,
                                 const float* szptr,
                                 const float* sxptr,
                                 const float* syptr,
                                 const float* rho,
                                 const float  dt,
                                 const float  dzi,
                                 const float  dxi,
                                 const float  dyi,
                                 const int    nz0,
                                 const int    nzf,
                                 const int    nx0,
                                 const int    nxf,
                                 const int    ny0,
                                 const int    nyf,
                                 const int    SZ,
                                 const int    SX,
                                 const int    SY,
                                 const int    dimmz,
                                 const int    dimmx,
                                 void*        stream);


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

/* CUDA STUB */
extern
void compute_component_vcell_BR_cuda ( float* vptr,
                                 const float* szptr,
                                 const float* sxptr,
                                 const float* syptr,
                                 const float* rho,
                                 const float  dt,
                                 const float  dzi,
                                 const float  dxi,
                                 const float  dyi,
                                 const int    nz0,
                                 const int    nzf,
                                 const int    nx0,
                                 const int    nxf,
                                 const int    ny0,
                                 const int    nyf,
                                 const int    SZ,
                                 const int    SX,
                                 const int    SY,
                                 const int    dimmz,
                                 const int    dimmx,
                                 void*        stream);


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

/* CUDA STUB */
extern
void compute_component_vcell_BL_cuda ( float* vptr,
                                 const float* szptr,
                                 const float* sxptr,
                                 const float* syptr,
                                 const float* rho,
                                 const float  dt,
                                 const float  dzi,
                                 const float  dxi,
                                 const float  dyi,
                                 const int    nz0,
                                 const int    nzf,
                                 const int    nx0,
                                 const int    nxf,
                                 const int    ny0,
                                 const int    nyf,
                                 const int    SZ,
                                 const int    SX,
                                 const int    SY,
                                 const int    dimmz,
                                 const int    dimmx,
                                 void*        stream);


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

#if defined(_OPENACC)
#pragma acc routine seq
#endif
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

#if defined(_OPENACC)
#pragma acc routine seq
#endif
real cell_coeff_BR ( const real* restrict ptr, 
                     const integer z, 
                     const integer x, 
                     const integer y, 
                     const integer dimmz, 
                     const integer dimmx );
#if defined(_OPENACC)
#pragma acc routine seq
#endif
real cell_coeff_TL ( const real* restrict ptr, 
                     const integer z, 
                     const integer x, 
                     const integer y, 
                     const integer dimmz, 
                     const integer dimmx );
#if defined(_OPENACC)
#pragma acc routine seq
#endif
real cell_coeff_BL ( const real* restrict ptr, 
                     const integer z, 
                     const integer x, 
                     const integer y, 
                     const integer dimmz, 
                     const integer dimmx );
#if defined(_OPENACC)
#pragma acc routine seq
#endif
real cell_coeff_TR ( const real* restrict ptr, 
                     const integer z, 
                     const integer x, 
                     const integer y, 
                     const integer dimmz, 
                     const integer dimmx );

#if defined(_OPENACC)
#pragma acc routine seq
#endif
real cell_coeff_ARTM_BR ( const real* restrict ptr, 
                          const integer z, 
                          const integer x, 
                          const integer y, 
                          const integer dimmz, 
                          const integer dimmx);
#if defined(_OPENACC)
#pragma acc routine seq
#endif
real cell_coeff_ARTM_TL ( const real* restrict ptr, 
                          const integer z, 
                          const integer x, 
                          const integer y, 
                          const integer dimmz, 
                          const integer dimmx);
#if defined(_OPENACC)
#pragma acc routine seq
#endif
real cell_coeff_ARTM_BL ( const real* restrict ptr, 
                          const integer z, 
                          const integer x, 
                          const integer y, 
                          const integer dimmz, 
                          const integer dimmx);
#if defined(_OPENACC)
#pragma acc routine seq
#endif
real cell_coeff_ARTM_TR ( const real* restrict ptr, 
                          const integer z, 
                          const integer x, 
                          const integer y, 
                          const integer dimmz, 
                          const integer dimmx);

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
                                 const phase_t phase);
/* STUB CUDA */
extern
void compute_component_scell_TR_cuda ( float* sxxptr,
                                       float* syyptr,
                                       float* szzptr,
                                       float* syzptr,
                                       float* sxzptr,
                                       float* sxyptr,
                                 const float* vxu,
                                 const float* vxv,
                                 const float* vxw,
                                 const float* vyu,
                                 const float* vyv,
                                 const float* vyw,
                                 const float* vzu,
                                 const float* vzv,
                                 const float* vzw,
                                 const float* cc11,
                                 const float* cc12,
                                 const float* cc13,
                                 const float* cc14,
                                 const float* cc15,
                                 const float* cc16,
                                 const float* cc22,
                                 const float* cc23,
                                 const float* cc24,
                                 const float* cc25,
                                 const float* cc26,
                                 const float* cc33,
                                 const float* cc34,
                                 const float* cc35,
                                 const float* cc36,
                                 const float* cc44,
                                 const float* cc45,
                                 const float* cc46,
                                 const float* cc55,
                                 const float* cc56,
                                 const float* cc66,
                                 const float  dt,
                                 const float  dzi,
                                 const float  dxi,
                                 const float  dyi,
                                 const int    nz0,
                                 const int    nzf,
                                 const int    nx0,
                                 const int    nxf,
                                 const int    ny0,
                                 const int    nyf,
                                 const int    SZ,
                                 const int    SX,
                                 const int    SY,
                                 const int    dimmz,
                                 const int    dimmx,
                                 void*        stream);


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
/* CUDA STUB */
extern
void compute_component_scell_TL_cuda ( float* sxxptr,
                                       float* syyptr,
                                       float* szzptr,
                                       float* syzptr,
                                       float* sxzptr,
                                       float* sxyptr,
                                 const float* vxu,
                                 const float* vxv,
                                 const float* vxw,
                                 const float* vyu,
                                 const float* vyv,
                                 const float* vyw,
                                 const float* vzu,
                                 const float* vzv,
                                 const float* vzw,
                                 const float* cc11,
                                 const float* cc12,
                                 const float* cc13,
                                 const float* cc14,
                                 const float* cc15,
                                 const float* cc16,
                                 const float* cc22,
                                 const float* cc23,
                                 const float* cc24,
                                 const float* cc25,
                                 const float* cc26,
                                 const float* cc33,
                                 const float* cc34,
                                 const float* cc35,
                                 const float* cc36,
                                 const float* cc44,
                                 const float* cc45,
                                 const float* cc46,
                                 const float* cc55,
                                 const float* cc56,
                                 const float* cc66,
                                 const float  dt,
                                 const float  dzi,
                                 const float  dxi,
                                 const float  dyi,
                                 const int    nz0,
                                 const int    nzf,
                                 const int    nx0,
                                 const int    nxf,
                                 const int    ny0,
                                 const int    nyf,
                                 const int    SZ,
                                 const int    SX,
                                 const int    SY,
                                 const int    dimmz,
                                 const int    dimmx,
                                 void*        stream);


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

/* CUDA STUB */
extern 
void compute_component_scell_BR_cuda ( float* sxxptr,
                                       float* syyptr,
                                       float* szzptr,
                                       float* syzptr,
                                       float* sxzptr,
                                       float* sxyptr,
                                 const float* vxu,
                                 const float* vxv,
                                 const float* vxw,
                                 const float* vyu,
                                 const float* vyv,
                                 const float* vyw,
                                 const float* vzu,
                                 const float* vzv,
                                 const float* vzw,
                                 const float* cc11,
                                 const float* cc12,
                                 const float* cc13,
                                 const float* cc14,
                                 const float* cc15,
                                 const float* cc16,
                                 const float* cc22,
                                 const float* cc23,
                                 const float* cc24,
                                 const float* cc25,
                                 const float* cc26,
                                 const float* cc33,
                                 const float* cc34,
                                 const float* cc35,
                                 const float* cc36,
                                 const float* cc44,
                                 const float* cc45,
                                 const float* cc46,
                                 const float* cc55,
                                 const float* cc56,
                                 const float* cc66,
                                 const float  dt,
                                 const float  dzi,
                                 const float  dxi,
                                 const float  dyi,
                                 const int    nz0,
                                 const int    nzf,
                                 const int    nx0,
                                 const int    nxf,
                                 const int    ny0,
                                 const int    nyf,
                                 const int    SZ,
                                 const int    SX,
                                 const int    SY,
                                 const int    dimmz,
                                 const int    dimmx,
                                 void*        stream);


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
/* CUDA STUB */
extern
void compute_component_scell_BL_cuda ( float* sxxptr,
                                       float* syyptr,
                                       float* szzptr,
                                       float* syzptr,
                                       float* sxzptr,
                                       float* sxyptr,
                                 const float* vxu,
                                 const float* vxv,
                                 const float* vxw,
                                 const float* vyu,
                                 const float* vyv,
                                 const float* vyw,
                                 const float* vzu,
                                 const float* vzv,
                                 const float* vzw,
                                 const float* cc11,
                                 const float* cc12,
                                 const float* cc13,
                                 const float* cc14,
                                 const float* cc15,
                                 const float* cc16,
                                 const float* cc22,
                                 const float* cc23,
                                 const float* cc24,
                                 const float* cc25,
                                 const float* cc26,
                                 const float* cc33,
                                 const float* cc34,
                                 const float* cc35,
                                 const float* cc36,
                                 const float* cc44,
                                 const float* cc45,
                                 const float* cc46,
                                 const float* cc55,
                                 const float* cc56,
                                 const float* cc66,
                                 const float  dt,
                                 const float  dzi,
                                 const float  dxi,
                                 const float  dyi,
                                 const int    nz0,
                                 const int    nzf,
                                 const int    nx0,
                                 const int    nxf,
                                 const int    ny0,
                                 const int    nyf,
                                 const int    SZ,
                                 const int    SX,
                                 const int    SY,
                                 const int    dimmz,
                                 const int    dimmx,
                                 void*        stream);


#ifdef __cplusplus
}
#endif /* extern "C" */

#endif /* end of _FWI_PROPAGATOR_H_ definition */
