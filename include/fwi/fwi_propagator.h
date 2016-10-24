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

typedef enum {back_offset, forw_offset} offset_t;
typedef enum {ONE_R, ONE_L, TWO, H2D, D2H} phase_t;

/*
 * \brief linearizes 3D indexing
 *
 * \param z     [in] index in Z dimension (contiguous in memory)
 * \param x     [in] index in X dimension (second most contiguous)
 * \param y     [in] index in Y dimension (less contiguous)
 * \param dimmz [in] size of Z dimension
 * \param dimmx [in] size of X dimension
 * \return linear index
 */
integer IDX (const integer z,
             const integer x,
             const integer y,
             const integer dimmz,
             const integer dimmx);

/*
 * \brief computes a 8-point stencil over Z dimension
 *
 * \param off   [in] offset (0 if BACKWARD or 1 if FORWARD)
 * \param ptr   [in] base pointer of a 3D volume
 * \param dzi   [in] delta in Z dimension
 * \param z     [in] index in Z dimension (contiguous in memory)
 * \param x     [in] index in X dimension (second most contiguous)
 * \param y     [in] index in Y dimension (less contiguous)
 * \param dimmz [in] size of Z dimension
 * \param dimmx [in] size of X dimension
 * \return stencil computation in Z dim
 */
real stencil_Z(const integer off,
               const real* restrict ptr,
               const real    dzi,
               const integer z,
               const integer x,
               const integer y,
               const integer dimmz,
               const integer dimmx);

/*
 * \brief computes a 8-point stencil over X dimension
 *
 * \param off   [in] offset (0 if BACKWARD or 1 if FORWARD)
 * \param ptr   [in] base pointer of a 3D volume
 * \param dxi   [in] delta in X dimension
 * \param z     [in] index in Z dimension (contiguous in memory)
 * \param x     [in] index in X dimension (second most contiguous)
 * \param y     [in] index in Y dimension (less contiguous)
 * \param dimmz [in] size of Z dimension
 * \param dimmx [in] size of X dimension
 * \return stencil computation in X dim
 */
real stencil_X(const integer off,
               const real* restrict ptr,
               const real dxi,
               const integer z,
               const integer x,
               const integer y,
               const integer dimmz,
               const integer dimmx);

/*
 * \brief computes a 8-point stencil over Y dimension
 *
 * \param off   [in] offset (0 if BACKWARD or 1 if FORWARD)
 * \param ptr   [in] base pointer of a 3D volume
 * \param dyi   [in] delta in Y dimension
 * \param z     [in] index in Z dimension (contiguous in memory)
 * \param x     [in] index in X dimension (second most contiguous)
 * \param y     [in] index in Y dimension (less contiguous)
 * \param dimmz [in] size of Z dimension
 * \param dimmx [in] size of X dimension
 * \return stencil computation in Y dim
 */
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
/*                             VELOCITY COMPUTATIONS                              */
/*                                                                                */
/* ------------------------------------------------------------------------------ */

/*
 * \brief computes rho BL coefficient
 *
 * \param rho   [in] base pointer of a 3D volume
 * \param z     [in] index in Z dimension (contiguous in memory)
 * \param x     [in] index in X dimension (second most contiguous)
 * \param y     [in] index in Y dimension (less contiguous)
 * \param dimmz [in] size of Z dimension
 * \param dimmx [in] size of X dimension
 * \return rho BL coefficient
 */
real rho_BL ( const real* restrict rho,
              const integer z,
              const integer x,
              const integer y,
              const integer dimmz,
              const integer dimmx);

/*
 * \brief computes rho TR coefficient
 *
 * \param rho   [in] base pointer of a 3D volume
 * \param z     [in] index in Z dimension (contiguous in memory)
 * \param x     [in] index in X dimension (second most contiguous)
 * \param y     [in] index in Y dimension (less contiguous)
 * \param dimmz [in] size of Z dimension
 * \param dimmx [in] size of X dimension
 * \return rho TR coefficient
 */
real rho_TR ( const real* restrict rho,
              const integer z,
              const integer x,
              const integer y,
              const integer dimmz,
              const integer dimmx);

/*
 * \brief computes rho BR coefficient
 *
 * \param rho   [in] base pointer of a 3D volume
 * \param z     [in] index in Z dimension (contiguous in memory)
 * \param x     [in] index in X dimension (second most contiguous)
 * \param y     [in] index in Y dimension (less contiguous)
 * \param dimmz [in] size of Z dimension
 * \param dimmx [in] size of X dimension
 * \return rho BR coefficient
 */
real rho_BR ( const real* restrict rho,
              const integer z,
              const integer x,
              const integer y,
              const integer dimmz,
              const integer dimmx);

/*
 * \brief computes rho TL coefficient
 *
 * \param rho   [in] base pointer of a 3D volume
 * \param z     [in] index in Z dimension (contiguous in memory)
 * \param x     [in] index in X dimension (second most contiguous)
 * \param y     [in] index in Y dimension (less contiguous)
 * \param dimmz [in] size of Z dimension
 * \param dimmx [in] size of X dimension
 * \return rho TL coefficient
 */
real rho_TL ( const real* restrict rho,
              const integer z,
              const integer x,
              const integer y,
              const integer dimmz,
              const integer dimmx);

/*
 * \brief computes velocity TL
 *
 * \param vptr  [inout] pointer from point_v_t
 * \param szptr [in] pointer from point_s_t
 * \param sxptr [in] pointer from point_s_t
 * \param syptr [in] pointer from point_s_t
 * \param rho   [in] pointer
 * \param dt    [in] delta time
 * \param dzi   [in] delta Z dimension
 * \param dxi   [in] delta X dimension
 * \param dyi   [in] delta Y dimension
 * \param nz0   [in] first element in Z dimension (range [nz0:nzf) )
 * \param nzf   [in] last  element in Z dimension
 * \param nx0   [in] first element in X dimension (range [nx0:nxf) )
 * \param nxf   [in] last  element in X dimension
 * \param ny0   [in] first element in Y dimension (range [ny0:nyf) )
 * \þaram nyf   [in] last  element in Y dimension
 * \param SZ    [in] offset in Z dimension
 * \param SX    [in] offset in X dimension
 * \param SY    [in] offset in Y dimension
 * \param dimmz [in] size of Z dimension
 * \param dimmx [in] size of X dimension
 * \param phase [in] enum that identifies the phase (useful for asynchronous exection with openacc)
 */
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
                                 const offset_t       SZ,
                                 const offset_t       SX,
                                 const offset_t       SY,
                                 const integer        dimmz,
                                 const integer        dimmx,
                                 const phase_t        phase);

#if defined(USE_CUDA)
//// COMPLETE CUDA STUB ////

/*
 * \brief CUDA implementation of 'compute_component_vcell_TL'
 *
 * This function acts as a 'glue' between 'compute_component_vcell_TL' and the CUDA Kernel
 *
 */
extern
void compute_component_vcell_TL_cuda ( /* COMPLETE ME */ );
////////////////////////////
#endif /* end USE_CUDA */

/*
 * \brief computes velocity TR
 *
 * \param vptr  [inout] pointer from point_v_t
 * \param szptr [in] pointer from point_s_t
 * \param sxptr [in] pointer from point_s_t
 * \param syptr [in] pointer from point_s_t
 * \param rho   [in] pointer
 * \param dt    [in] delta time
 * \param dzi   [in] delta Z dimension
 * \param dxi   [in] delta X dimension
 * \param dyi   [in] delta Y dimension
 * \param nz0   [in] first element in Z dimension (range [nz0:nzf) )
 * \param nzf   [in] last  element in Z dimension
 * \param nx0   [in] first element in X dimension (range [nx0:nxf) )
 * \param nxf   [in] last  element in X dimension
 * \param ny0   [in] first element in Y dimension (range [ny0:nyf) )
 * \þaram nyf   [in] last  element in Y dimension
 * \param SZ    [in] offset in Z dimension
 * \param SX    [in] offset in X dimension
 * \param SY    [in] offset in Y dimension
 * \param dimmz [in] size of Z dimension
 * \param dimmx [in] size of X dimension
 * \param phase [in] enum that identifies the phase (useful for asynchronous exection with openacc)
 */
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
                                 const offset_t       SZ,
                                 const offset_t       SX,
                                 const offset_t       SY,
                                 const integer        dimmz,
                                 const integer        dimmx,
                                 const phase_t        phase);

#if defined(USE_CUDA)
//// COMPLETE CUDA STUB ////

/*
 * \brief CUDA implementation of 'compute_component_vcell_TR'
 *
 * This function acts as a 'glue' between 'compute_component_vcell_TR' and the CUDA Kernel
 *
 */
extern
void compute_component_vcell_TR_cuda ( /* COMPLETE ME */ );
////////////////////////////
#endif /* end USE_CUDA */

/*
 * \brief computes velocity BR
 *
 * \param vptr  [inout] pointer from point_v_t
 * \param szptr [in] pointer from point_s_t
 * \param sxptr [in] pointer from point_s_t
 * \param syptr [in] pointer from point_s_t
 * \param rho   [in] pointer
 * \param dt    [in] delta time
 * \param dzi   [in] delta Z dimension
 * \param dxi   [in] delta X dimension
 * \param dyi   [in] delta Y dimension
 * \param nz0   [in] first element in Z dimension (range [nz0:nzf) )
 * \param nzf   [in] last  element in Z dimension
 * \param nx0   [in] first element in X dimension (range [nx0:nxf) )
 * \param nxf   [in] last  element in X dimension
 * \param ny0   [in] first element in Y dimension (range [ny0:nyf) )
 * \þaram nyf   [in] last  element in Y dimension
 * \param SZ    [in] offset in Z dimension
 * \param SX    [in] offset in X dimension
 * \param SY    [in] offset in Y dimension
 * \param dimmz [in] size of Z dimension
 * \param dimmx [in] size of X dimension
 * \param phase [in] enum that identifies the phase (useful for asynchronous exection with openacc)
 */
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
                                 const offset_t       SZ,
                                 const offset_t       SX,
                                 const offset_t       SY,
                                 const integer        dimmz,
                                 const integer        dimmx,
                                 const phase_t        phase);

#if defined(USE_CUDA)
//// COMPLETE CUDA STUB ////

/*
 * \brief CUDA implementation of 'compute_component_vcell_BR'
 *
 * This function acts as a 'glue' between 'compute_component_vcell_BR' and the CUDA Kernel
 *
 */
extern
void compute_component_vcell_BR_cuda ( /* COMPLETE ME */ );
////////////////////////////
#endif /* end USE_CUDA */

/*
 * \brief computes velocity BL
 *
 * \param vptr  [inout] pointer from point_v_t
 * \param szptr [in] pointer from point_s_t
 * \param sxptr [in] pointer from point_s_t
 * \param syptr [in] pointer from point_s_t
 * \param rho   [in] pointer
 * \param dt    [in] delta time
 * \param dzi   [in] delta Z dimension
 * \param dxi   [in] delta X dimension
 * \param dyi   [in] delta Y dimension
 * \param nz0   [in] first element in Z dimension (range [nz0:nzf) )
 * \param nzf   [in] last  element in Z dimension
 * \param nx0   [in] first element in X dimension (range [nx0:nxf) )
 * \param nxf   [in] last  element in X dimension
 * \param ny0   [in] first element in Y dimension (range [ny0:nyf) )
 * \þaram nyf   [in] last  element in Y dimension
 * \param SZ    [in] offset in Z dimension
 * \param SX    [in] offset in X dimension
 * \param SY    [in] offset in Y dimension
 * \param dimmz [in] size of Z dimension
 * \param dimmx [in] size of X dimension
 * \param phase [in] enum that identifies the phase (useful for asynchronous exection with openacc)
 */
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
                                 const offset_t       SZ,
                                 const offset_t       SX,
                                 const offset_t       SY,
                                 const integer        dimmz,
                                 const integer        dimmx,
                                 const phase_t        phase);

#if defined(USE_CUDA)
//// COMPLETE CUDA STUB ////

/*
 * \brief CUDA implementation of 'compute_component_vcell_BL'
 *
 * This function acts as a 'glue' between 'compute_component_vcell_BL' and the CUDA Kernel
 *
 */
extern
void compute_component_vcell_BL_cuda ( /* COMPLETE ME */ );
////////////////////////////
#endif /* end USE_CUDA */


/*
 * \brief computes velocity propagation
 *
 * \param v     [in] v_t structure
 * \param s     [in] s_t structure
 * \param coeffs [in] coeff_t coefficients (not used)
 * \param rho   [in] pointer
 * \param dt    [in] delta time
 * \param dzi   [in] delta Z dimension
 * \param dxi   [in] delta X dimension
 * \param dyi   [in] delta Y dimension
 * \param nz0   [in] first element in Z dimension (range [nz0:nzf) )
 * \param nzf   [in] last  element in Z dimension
 * \param nx0   [in] first element in X dimension (range [nx0:nxf) )
 * \param nxf   [in] last  element in X dimension
 * \param ny0   [in] first element in Y dimension (range [ny0:nyf) )
 * \þaram nyf   [in] last  element in Y dimension
 * \param dimmz [in] size of Z dimension
 * \param dimmx [in] size of X dimension
 * \param phase [in] enum that identifies the phase (useful for asynchronous exection with openacc)
 */
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
/*                              STRESS COMPUTATIONS                               */
/*                                                                                */
/* ------------------------------------------------------------------------------ */

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

real cell_coeff_BR ( const real* restrict ptr,
                     const integer z,
                     const integer x,
                     const integer y,
                     const integer dimmz,
                     const integer dimmx );

real cell_coeff_TL ( const real* restrict ptr,
                     const integer z,
                     const integer x,
                     const integer y,
                     const integer dimmz,
                     const integer dimmx );

real cell_coeff_BL ( const real* restrict ptr,
                     const integer z,
                     const integer x,
                     const integer y,
                     const integer dimmz,
                     const integer dimmx );

real cell_coeff_TR ( const real* restrict ptr,
                     const integer z,
                     const integer x,
                     const integer y,
                     const integer dimmz,
                     const integer dimmx );

real cell_coeff_ARTM_BR ( const real* restrict ptr,
                          const integer z,
                          const integer x,
                          const integer y,
                          const integer dimmz,
                          const integer dimmx);

real cell_coeff_ARTM_TL ( const real* restrict ptr,
                          const integer z,
                          const integer x,
                          const integer y,
                          const integer dimmz,
                          const integer dimmx);

real cell_coeff_ARTM_BL ( const real* restrict ptr,
                          const integer z,
                          const integer x,
                          const integer y,
                          const integer dimmz,
                          const integer dimmx);

real cell_coeff_ARTM_TR ( const real* restrict ptr,
                          const integer z,
                          const integer x,
                          const integer y,
                          const integer dimmz,
                          const integer dimmx);

/*
 * \brief computes stresses TR
 *
 * \param v       [in] s_t structure
 * \param vnode_z [in] pointer from point_v_t
 * \param vnode_x [in] pointer from point_v_t
 * \param vnode_y [in] pointer from point_v_t
 * \param dt      [in] delta time
 * \param dzi     [in] delta Z dimension
 * \param dxi     [in] delta X dimension
 * \param dyi     [in] delta Y dimension
 * \param nz0     [in] first element in Z dimension (range [nz0:nzf) )
 * \param nzf     [in] last  element in Z dimension
 * \param nx0     [in] first element in X dimension (range [nx0:nxf) )
 * \param nxf     [in] last  element in X dimension
 * \param ny0     [in] first element in Y dimension (range [ny0:nyf) )
 * \þaram nyf     [in] last  element in Y dimension
 * \param SZ      [in] offset in Z dimension
 * \param SX      [in] offset in X dimension
 * \param SY      [in] offset in Y dimension
 * \param dimmz   [in] size of Z dimension
 * \param dimmx   [in] size of X dimension
 * \param phase   [in] enum that identifies the phase (useful for asynchronous exection with openacc)
 */
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
                                 const integer  dimmx,
                                 const phase_t  phase);

#if defined(USE_CUDA)
//// COMPLETE CUDA STUB ////

/*
 * \brief CUDA implementation of 'compute_component_scell_TR'
 *
 * This function acts as a 'glue' between 'compute_component_scell_TR' and the CUDA Kernel
 *
 */
extern
void compute_component_scell_TR_cuda ( /* COMPLETE ME */ );
////////////////////////////
#endif /* end USE_CUDA */

/*
 * \brief computes stresses TL
 *
 * \param v       [in] s_t structure
 * \param vnode_z [in] pointer from point_v_t
 * \param vnode_x [in] pointer from point_v_t
 * \param vnode_y [in] pointer from point_v_t
 * \param dt      [in] delta time
 * \param dzi     [in] delta Z dimension
 * \param dxi     [in] delta X dimension
 * \param dyi     [in] delta Y dimension
 * \param nz0     [in] first element in Z dimension (range [nz0:nzf) )
 * \param nzf     [in] last  element in Z dimension
 * \param nx0     [in] first element in X dimension (range [nx0:nxf) )
 * \param nxf     [in] last  element in X dimension
 * \param ny0     [in] first element in Y dimension (range [ny0:nyf) )
 * \þaram nyf     [in] last  element in Y dimension
 * \param SZ      [in] offset in Z dimension
 * \param SX      [in] offset in X dimension
 * \param SY      [in] offset in Y dimension
 * \param dimmz   [in] size of Z dimension
 * \param dimmx   [in] size of X dimension
 * \param phase   [in] enum that identifies the phase (useful for asynchronous exection with openacc)
 */
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
                                 const integer  dimmx,
                                 const phase_t  phase);

#if defined(USE_CUDA)
//// COMPLETE CUDA STUB ////

/*
 * \brief CUDA implementation of 'compute_component_scell_TL'
 *
 * This function acts as a 'glue' between 'compute_component_scell_TL' and the CUDA Kernel
 *
 */
extern
void compute_component_scell_TL_cuda ( /* COMPLETE ME */ );
////////////////////////////
#endif /* end USE_CUDA */

/*
 * \brief computes stresses BR
 *
 * \param v       [in] s_t structure
 * \param vnode_z [in] pointer from point_v_t
 * \param vnode_x [in] pointer from point_v_t
 * \param vnode_y [in] pointer from point_v_t
 * \param dt      [in] delta time
 * \param dzi     [in] delta Z dimension
 * \param dxi     [in] delta X dimension
 * \param dyi     [in] delta Y dimension
 * \param nz0     [in] first element in Z dimension (range [nz0:nzf) )
 * \param nzf     [in] last  element in Z dimension
 * \param nx0     [in] first element in X dimension (range [nx0:nxf) )
 * \param nxf     [in] last  element in X dimension
 * \param ny0     [in] first element in Y dimension (range [ny0:nyf) )
 * \þaram nyf     [in] last  element in Y dimension
 * \param SZ      [in] offset in Z dimension
 * \param SX      [in] offset in X dimension
 * \param SY      [in] offset in Y dimension
 * \param dimmz   [in] size of Z dimension
 * \param dimmx   [in] size of X dimension
 * \param phase   [in] enum that identifies the phase (useful for asynchronous exection with openacc)
 */
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
                                 const integer  dimmx,
                                 const phase_t  phase);
#if defined(USE_CUDA)
//// COMPLETE CUDA STUB ////

/*
 * \brief CUDA implementation of 'compute_component_scell_BR'
 *
 * This function acts as a 'glue' between 'compute_component_scell_BR' and the CUDA Kernel
 *
 */
extern
void compute_component_scell_BR_cuda ( /* COMPLETE ME */ );
////////////////////////////
#endif /* end USE_CUDA */

/*
 * \brief computes stresses BL
 *
 * \param v       [in] s_t structure
 * \param vnode_z [in] pointer from point_v_t
 * \param vnode_x [in] pointer from point_v_t
 * \param vnode_y [in] pointer from point_v_t
 * \param dt      [in] delta time
 * \param dzi     [in] delta Z dimension
 * \param dxi     [in] delta X dimension
 * \param dyi     [in] delta Y dimension
 * \param nz0     [in] first element in Z dimension (range [nz0:nzf) )
 * \param nzf     [in] last  element in Z dimension
 * \param nx0     [in] first element in X dimension (range [nx0:nxf) )
 * \param nxf     [in] last  element in X dimension
 * \param ny0     [in] first element in Y dimension (range [ny0:nyf) )
 * \þaram nyf     [in] last  element in Y dimension
 * \param SZ      [in] offset in Z dimension
 * \param SX      [in] offset in X dimension
 * \param SY      [in] offset in Y dimension
 * \param dimmz   [in] size of Z dimension
 * \param dimmx   [in] size of X dimension
 * \param phase   [in] enum that identifies the phase (useful for asynchronous exection with openacc)
 */
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
                                 const integer  dimmx,
                                 const phase_t  phase);
#if defined(USE_CUDA)
//// COMPLETE CUDA STUB ////

/*
 * \brief CUDA implementation of 'compute_component_scell_BL'
 *
 * This function acts as a 'glue' between 'compute_component_scell_BL' and the CUDA Kernel
 *
 */
extern
void compute_component_scell_BL_cuda ( /* COMPLETE ME */ );
////////////////////////////
#endif /* end USE_CUDA */


#ifdef __cplusplus
}
#endif /* extern "C" */

#endif /* end of _FWI_PROPAGATOR_H_ definition */
