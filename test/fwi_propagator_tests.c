#include <unity.h>
#include <unity_fixture.h>

#include "fwi_kernel.h"
#include "fwi_propagator.h"

TEST_GROUP(propagator);

// Helper functions
void init_array( real* restrict array, const integer length );
void copy_array( real* restrict dest, real* restrict src, const integer length );


////// SetUp/SetDown /////
const integer dimmz = 32;
const integer dimmx = 16;
const integer dimmy = 16;
      integer nelems;

v_t v_ref;
s_t s_ref;
coeff_t c_ref;
real* rho_ref;

v_t v_cal;
s_t s_cal;
coeff_t c_cal;
real* rho_cal;

void init_array( real* restrict array, const integer length )
{
    for (integer i = 0; i < length; i++)
        array[i] = (i/(1.0+i)) + (rand() + 1.0) / (1.0 * RAND_MAX);
}

void copy_array( real* restrict dest, real* restrict src, const integer length )
{
    for (integer i = 0; i < length; i++)
        dest[i] = src[i];
}


TEST_SETUP(propagator)
{
    nelems = dimmz * dimmx * dimmy;

    alloc_memory_shot(nelems, &c_ref, &s_ref, &v_ref, &rho_ref);
    alloc_memory_shot(nelems, &c_cal, &s_cal, &v_cal, &rho_cal);

    /* constants initialization */
    init_array(c_ref.c11, nelems);
    init_array(c_ref.c12, nelems);
    init_array(c_ref.c13, nelems);
    init_array(c_ref.c14, nelems);
    init_array(c_ref.c15, nelems);
    init_array(c_ref.c16, nelems);
    
    init_array(c_ref.c22, nelems);
    init_array(c_ref.c23, nelems);
    init_array(c_ref.c24, nelems);
    init_array(c_ref.c25, nelems);
    init_array(c_ref.c26, nelems);

    init_array(c_ref.c33, nelems);
    init_array(c_ref.c34, nelems);
    init_array(c_ref.c35, nelems);
    init_array(c_ref.c36, nelems);

    init_array(c_ref.c44, nelems);
    init_array(c_ref.c45, nelems);
    init_array(c_ref.c46, nelems);

    init_array(c_ref.c55, nelems);
    init_array(c_ref.c56, nelems);

    init_array(c_ref.c66, nelems);

    /* velocity [reference] initialization */
    init_array(v_ref.tl.u, nelems);
    init_array(v_ref.tl.v, nelems);
    init_array(v_ref.tl.w, nelems);

    init_array(v_ref.tr.u, nelems);
    init_array(v_ref.tr.v, nelems);
    init_array(v_ref.tr.w, nelems);
    
    init_array(v_ref.bl.u, nelems);
    init_array(v_ref.bl.v, nelems);
    init_array(v_ref.bl.w, nelems);

    init_array(v_ref.br.u, nelems);
    init_array(v_ref.br.v, nelems);
    init_array(v_ref.br.w, nelems);

    /* stresses [reference] initialization */
    init_array(s_ref.tl.zz, nelems);
    init_array(s_ref.tl.xz, nelems);
    init_array(s_ref.tl.yz, nelems);
    init_array(s_ref.tl.xx, nelems);
    init_array(s_ref.tl.xy, nelems);
    init_array(s_ref.tl.yy, nelems);

    init_array(s_ref.tr.zz, nelems);
    init_array(s_ref.tr.xz, nelems);
    init_array(s_ref.tr.yz, nelems);
    init_array(s_ref.tr.xx, nelems);
    init_array(s_ref.tr.xy, nelems);
    init_array(s_ref.tr.yy, nelems);

    init_array(s_ref.bl.zz, nelems);
    init_array(s_ref.bl.xz, nelems);
    init_array(s_ref.bl.yz, nelems);
    init_array(s_ref.bl.xx, nelems);
    init_array(s_ref.bl.xy, nelems);
    init_array(s_ref.bl.yy, nelems);

    init_array(s_ref.br.zz, nelems);
    init_array(s_ref.br.xz, nelems);
    init_array(s_ref.br.yz, nelems);
    init_array(s_ref.br.xx, nelems);
    init_array(s_ref.br.xy, nelems);
    init_array(s_ref.br.yy, nelems);

    /* init rho */
    init_array(rho_ref, nelems);

    ///////////////////////////////
    // Copy 'ref' values into 'cal'
    // to start from equal solutions

    copy_array(v_cal.tl.u, v_ref.tl.u, nelems);
    copy_array(v_cal.tl.v, v_ref.tl.v, nelems);
    copy_array(v_cal.tl.w, v_ref.tl.w, nelems);
           
    copy_array(v_cal.tr.u, v_ref.tr.u, nelems);
    copy_array(v_cal.tr.v, v_ref.tr.v, nelems);
    copy_array(v_cal.tr.w, v_ref.tr.w, nelems);
    
    copy_array(v_cal.bl.u, v_ref.bl.u, nelems);
    copy_array(v_cal.bl.v, v_ref.bl.v, nelems);
    copy_array(v_cal.bl.w, v_ref.bl.w, nelems);
    
    copy_array(v_cal.br.u, v_ref.br.u, nelems);
    copy_array(v_cal.br.v, v_ref.br.v, nelems);
    copy_array(v_cal.br.w, v_ref.br.w, nelems);

    copy_array(s_cal.tl.zz, s_ref.tl.zz, nelems);
    copy_array(s_cal.tl.xz, s_ref.tl.xz, nelems);
    copy_array(s_cal.tl.yz, s_ref.tl.yz, nelems);
    copy_array(s_cal.tl.xx, s_ref.tl.xx, nelems);
    copy_array(s_cal.tl.xy, s_ref.tl.xy, nelems);
    copy_array(s_cal.tl.yy, s_ref.tl.yy, nelems);
     
    copy_array(s_cal.tr.zz, s_ref.tr.zz, nelems);
    copy_array(s_cal.tr.xz, s_ref.tr.xz, nelems);
    copy_array(s_cal.tr.yz, s_ref.tr.yz, nelems);
    copy_array(s_cal.tr.xx, s_ref.tr.xx, nelems);
    copy_array(s_cal.tr.xy, s_ref.tr.xy, nelems);
    copy_array(s_cal.tr.yy, s_ref.tr.yy, nelems);
    
    copy_array(s_cal.bl.zz, s_ref.bl.zz, nelems);
    copy_array(s_cal.bl.xz, s_ref.bl.xz, nelems);
    copy_array(s_cal.bl.yz, s_ref.bl.yz, nelems);
    copy_array(s_cal.bl.xx, s_ref.bl.xx, nelems);
    copy_array(s_cal.bl.xy, s_ref.bl.xy, nelems);
    copy_array(s_cal.bl.yy, s_ref.bl.yy, nelems);
    
    copy_array(s_cal.br.zz, s_ref.br.zz, nelems);
    copy_array(s_cal.br.xz, s_ref.br.xz, nelems);
    copy_array(s_cal.br.yz, s_ref.br.yz, nelems);
    copy_array(s_cal.br.xx, s_ref.br.xx, nelems);
    copy_array(s_cal.br.xy, s_ref.br.xy, nelems);
    copy_array(s_cal.br.yy, s_ref.br.yy, nelems);
}

TEST_TEAR_DOWN(propagator)
{
    free_memory_shot(&c_ref, &s_ref, &v_ref, &rho_ref);
    free_memory_shot(&c_cal, &s_cal, &v_cal, &rho_cal);
}

/////// TESTS ////////////

TEST(propagator, IDX)
{
    TEST_ASSERT_EQUAL_INT( 0*dimmx*dimmz + 0*dimmz + 0, IDX(0, 0, 0, dimmz, dimmx) );
    TEST_ASSERT_EQUAL_INT( 1*dimmx*dimmz + 0*dimmz + 0, IDX(0, 0, 1, dimmz, dimmx) );
    TEST_ASSERT_EQUAL_INT( 0*dimmx*dimmz + 1*dimmz + 0, IDX(0, 1, 0, dimmz, dimmx) );
    TEST_ASSERT_EQUAL_INT( 1*dimmx*dimmz + 1*dimmz + 0, IDX(0, 1, 1, dimmz, dimmx) );
    TEST_ASSERT_EQUAL_INT( 0*dimmx*dimmz + 0*dimmz + 1, IDX(1, 0, 0, dimmz, dimmx) );
    TEST_ASSERT_EQUAL_INT( 1*dimmx*dimmz + 0*dimmz + 1, IDX(1, 0, 1, dimmz, dimmx) );
    TEST_ASSERT_EQUAL_INT( 0*dimmx*dimmz + 1*dimmz + 1, IDX(1, 1, 0, dimmz, dimmx) );
    TEST_ASSERT_EQUAL_INT( 1*dimmx*dimmz + 1*dimmz + 1, IDX(1, 1, 1, dimmz, dimmx) );
}

TEST(propagator, stencil_Z)
{
    const integer FORWARD  = 1;
    const integer BACKWARD = 0;

    const real dzi = 1.0;

    for (integer y =    0; y < dimmy;        y++)
    for (integer x =    0; x < dimmx;        x++)
    for (integer z = HALO; z < dimmz - HALO; z++)
    {
        TEST_ASSERT_EQUAL_FLOAT(stencil_Z(FORWARD,  v_ref.tl.u, dzi, z,   x, y, dimmz, dimmx),
                                stencil_Z(BACKWARD, v_ref.tl.u, dzi, z+1, x, y, dimmz, dimmx));
    }
}

TEST(propagator, stencil_X)
{
    const integer FORWARD  = 1;
    const integer BACKWARD = 0;

    const real dzi = 1.0;

    for (integer y =    0; y < dimmy;        y++)
    for (integer x = HALO; x < dimmx - HALO; x++)
    for (integer z =    0; z < dimmz;        z++)
    {
        TEST_ASSERT_EQUAL_FLOAT(stencil_X(FORWARD,  v_ref.tl.v, dzi, z, x,   y, dimmz, dimmx),
                                stencil_X(BACKWARD, v_ref.tl.v, dzi, z, x+1, y, dimmz, dimmx));
    }
}

TEST(propagator, stencil_Y)
{
    const integer FORWARD  = 1;
    const integer BACKWARD = 0;

    const real dzi = 1.0;

    for (integer y = HALO; y < dimmy - HALO; y++)
    for (integer x =    0; x < dimmx;        x++)
    for (integer z =    0; z < dimmz;        z++)
    {
        TEST_ASSERT_EQUAL_FLOAT(stencil_Y(FORWARD,  v_ref.tl.w, dzi, z, x, y,   dimmz, dimmx),
                                stencil_Y(BACKWARD, v_ref.tl.w, dzi, z, x, y+1, dimmz, dimmx));
    }
}

TEST(propagator, rho_BL)
{
    for (integer y = 0; y < dimmy;   y++)
    for (integer x = 0; x < dimmx;   x++)
    for (integer z = 0; z < dimmz-1; z++)
    {
        TEST_ASSERT_EQUAL_FLOAT((2.0f/(v_ref.bl.u[IDX(z  ,x,y,dimmz,dimmx)]+v_ref.bl.u[IDX(z+1,x,y,dimmz,dimmx)])),
                                rho_BL(v_ref.bl.u, z, x, y, dimmz, dimmx));
    }
}

TEST(propagator, rho_TR)
{
    for (integer y = 0; y < dimmy;   y++)
    for (integer x = 0; x < dimmx-1; x++)
    for (integer z = 0; z < dimmz;   z++)
    {
        TEST_ASSERT_EQUAL_FLOAT((2.0f/(v_ref.tr.u[IDX(z,x,y,dimmz,dimmx)]+v_ref.tr.u[IDX(z,x+1,y,dimmz,dimmx)])),
                                rho_TR(v_ref.tr.u, z, x, y, dimmz, dimmx));
    }
}

TEST(propagator, rho_BR)
{
    for (integer y = 0; y < dimmy-1; y++)
    for (integer x = 0; x < dimmx-1; x++)
    for (integer z = 0; z < dimmz-1; z++)
    {
        const real ref = ( 8.0f / ( v_ref.br.u[IDX(z  ,x  ,y  ,dimmz,dimmx)] +
                                    v_ref.br.u[IDX(z+1,x  ,y  ,dimmz,dimmx)] +
                                    v_ref.br.u[IDX(z  ,x+1,y  ,dimmz,dimmx)] +
                                    v_ref.br.u[IDX(z  ,x  ,y+1,dimmz,dimmx)] +
                                    v_ref.br.u[IDX(z  ,x+1,y+1,dimmz,dimmx)] +
                                    v_ref.br.u[IDX(z+1,x+1,y  ,dimmz,dimmx)] +
                                    v_ref.br.u[IDX(z+1,x  ,y+1,dimmz,dimmx)] +
                                    v_ref.br.u[IDX(z+1,x+1,y+1,dimmz,dimmx)]) );
        
        const real cal = rho_BR(v_ref.br.u, z, x, y, dimmz, dimmx);

        TEST_ASSERT_EQUAL_FLOAT(ref, cal);
    }
}

TEST(propagator, rho_TL)
{
    for (integer y = 0; y < dimmy-1; y++)
    for (integer x = 0; x < dimmx;   x++)
    for (integer z = 0; z < dimmz;   z++)
    {
        TEST_ASSERT_EQUAL_FLOAT((2.0f/(v_ref.tl.u[IDX(z,x,y,dimmz,dimmx)]+v_ref.tl.u[IDX(z,x,y+1,dimmz,dimmx)])),
                                rho_TL(v_ref.tl.u, z, x, y, dimmz, dimmx));
    }
}



TEST(propagator, compute_component_vcell_TL)
{
    const real     dt  = 1.0;
    const real     dzi = 1.0;
    const real     dxi = 1.0;
    const real     dyi = 1.0;
    const integer  nz0 = HALO;
    const integer  nzf = dimmz-HALO;
    const integer  nx0 = HALO;
    const integer  nxf = dimmx-HALO;
    const integer  ny0 = HALO;
    const integer  nyf = dimmy-HALO;
    const offset_t SZ = 0;
    const offset_t SX = 0;
    const offset_t SY = 0;
    const phase_t  phase = TWO;

    const real*    szptr = s_ref.bl.xz;
    const real*    sxptr = s_ref.tr.xx;
    const real*    syptr = s_ref.tl.xy;

    // REFERENCE CALCULATION -DON'T TOUCH-
    for(integer y=ny0; y < nyf; y++)
    {
        for(integer x=nx0; x < nxf; x++)
        {
            for(integer z=nz0; z < nzf; z++)
            {
                const real lrho = rho_TL(rho_ref, z, x, y, dimmz, dimmx);
                
                const real stx  = stencil_X( SX, sxptr, dxi, z, x, y, dimmz, dimmx);
                const real sty  = stencil_Y( SY, syptr, dyi, z, x, y, dimmz, dimmx);
                const real stz  = stencil_Z( SZ, szptr, dzi, z, x, y, dimmz, dimmx);
                
                v_ref.tl.u[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
    ///////////////////////////////////////

#if defined(_OPENACC)
    const integer start  = ((nzf-nz0) + 2*HALO) * ((nxf-nx0) + 2*HALO) * (ny0 - HALO);
    const integer end    = ((nzf-nz0) + 2*HALO) * ((nxf-nx0) + 2*HALO) * (nyf + HALO);
    const integer nelems = end - start;

    #pragma acc data copy(calculated[start:nelems]) \
                     copyin(szptr[start:nelems], sxptr[start:nelems], syptr[start:nelems], rho[start:nelems])
    {
#endif
        compute_component_vcell_TL( v_cal.tl.u, szptr, sxptr, syptr, rho_ref, 
            dt, dzi, dxi, dyi, 
            nz0, nzf, nx0, nxf, ny0, nyf, 
            SZ, SX, SY, dimmz, dimmx, phase);
#if defined(_OPENACC)
    }
#endif

    TEST_ASSERT_EQUAL_FLOAT_ARRAY( v_ref.tl.u, v_cal.tl.u, nelems );
}

TEST(propagator, compute_component_vcell_TR)
{
    const real     dt  = 1.0;
    const real     dzi = 1.0;
    const real     dxi = 1.0;
    const real     dyi = 1.0;
    const integer  nz0 = HALO;
    const integer  nzf = dimmz-HALO;
    const integer  nx0 = HALO;
    const integer  nxf = dimmx-HALO;
    const integer  ny0 = HALO;
    const integer  nyf = dimmy-HALO;
    const offset_t SZ = 0;
    const offset_t SX = 0;
    const offset_t SY = 0;
    const phase_t  phase = TWO;

    const real*    szptr = s_ref.br.xz;
    const real*    sxptr = s_ref.tl.xx;
    const real*    syptr = s_ref.tr.xy;

    // REFERENCE CALCULATION -DON'T TOUCH-
    for(integer y=ny0; y < nyf; y++)
    {
        for(integer x=nx0; x < nxf; x++)
        {
            for(integer z=nz0; z < nzf; z++)
            {
                const real lrho = rho_TR(rho_ref, z, x, y, dimmz, dimmx);
                
                const real stx  = stencil_X( SX, sxptr, dxi, z, x, y, dimmz, dimmx);
                const real sty  = stencil_Y( SY, syptr, dyi, z, x, y, dimmz, dimmx);
                const real stz  = stencil_Z( SZ, szptr, dzi, z, x, y, dimmz, dimmx);
                
                v_ref.tr.u[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
    ///////////////////////////////////////

#if defined(_OPENACC)
    const integer start  = ((nzf-nz0) + 2*HALO) * ((nxf-nx0) + 2*HALO) * (ny0 - HALO);
    const integer end    = ((nzf-nz0) + 2*HALO) * ((nxf-nx0) + 2*HALO) * (nyf + HALO);
    const integer nelems = end - start;

    #pragma acc data copy(calculated[start:nelems]) \
                     copyin(szptr[start:nelems], sxptr[start:nelems], syptr[start:nelems], rho[start:nelems])
    {
#endif
        compute_component_vcell_TR( v_cal.tr.u, szptr, sxptr, syptr, rho_ref, 
            dt, dzi, dxi, dyi, 
            nz0, nzf, nx0, nxf, ny0, nyf, 
            SZ, SX, SY, dimmz, dimmx, phase);
#if defined(_OPENACC)
    }
#endif

    TEST_ASSERT_EQUAL_FLOAT_ARRAY( v_ref.tr.u, v_cal.tr.u, nelems );
}


TEST(propagator, compute_component_vcell_BR)
{
    const real     dt  = 1.0;
    const real     dzi = 1.0;
    const real     dxi = 1.0;
    const real     dyi = 1.0;
    const integer  nz0 = HALO;
    const integer  nzf = dimmz-HALO;
    const integer  nx0 = HALO;
    const integer  nxf = dimmx-HALO;
    const integer  ny0 = HALO;
    const integer  nyf = dimmy-HALO;
    const offset_t SZ = 0;
    const offset_t SX = 0;
    const offset_t SY = 0;
    const phase_t  phase = TWO;

    const real*    szptr = s_ref.tr.xz;
    const real*    sxptr = s_ref.bl.xx;
    const real*    syptr = s_ref.br.xy;

    // REFERENCE CALCULATION -DON'T TOUCH-
    for(integer y=ny0; y < nyf; y++)
    {
        for(integer x=nx0; x < nxf; x++)
        {
            for(integer z=nz0; z < nzf; z++)
            {
                const real lrho = rho_BR(rho_ref, z, x, y, dimmz, dimmx);
                
                const real stx  = stencil_X( SX, sxptr, dxi, z, x, y, dimmz, dimmx);
                const real sty  = stencil_Y( SY, syptr, dyi, z, x, y, dimmz, dimmx);
                const real stz  = stencil_Z( SZ, szptr, dzi, z, x, y, dimmz, dimmx);
                
                v_ref.br.u[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
    ///////////////////////////////////////

#if defined(_OPENACC)
    const integer start  = ((nzf-nz0) + 2*HALO) * ((nxf-nx0) + 2*HALO) * (ny0 - HALO);
    const integer end    = ((nzf-nz0) + 2*HALO) * ((nxf-nx0) + 2*HALO) * (nyf + HALO);
    const integer nelems = end - start;

    #pragma acc data copy(calculated[start:nelems]) \
                     copyin(szptr[start:nelems], sxptr[start:nelems], syptr[start:nelems], rho[start:nelems])
    {
#endif
        compute_component_vcell_BR( v_cal.br.u, szptr, sxptr, syptr, rho_ref, 
            dt, dzi, dxi, dyi, 
            nz0, nzf, nx0, nxf, ny0, nyf, 
            SZ, SX, SY, dimmz, dimmx, phase);
#if defined(_OPENACC)
    }
#endif

    TEST_ASSERT_EQUAL_FLOAT_ARRAY( v_ref.br.u, v_cal.br.u, nelems );
}

TEST(propagator, compute_component_vcell_BL)
{
    const real     dt  = 1.0;
    const real     dzi = 1.0;
    const real     dxi = 1.0;
    const real     dyi = 1.0;
    const integer  nz0 = HALO;
    const integer  nzf = dimmz-HALO;
    const integer  nx0 = HALO;
    const integer  nxf = dimmx-HALO;
    const integer  ny0 = HALO;
    const integer  nyf = dimmy-HALO;
    const offset_t SZ = 0;
    const offset_t SX = 0;
    const offset_t SY = 0;
    const phase_t  phase = TWO;

    const real*    szptr = s_ref.tl.xz;
    const real*    sxptr = s_ref.br.xx;
    const real*    syptr = s_ref.bl.xy;

    // REFERENCE CALCULATION -DON'T TOUCH-
    for(integer y=ny0; y < nyf; y++)
    {
        for(integer x=nx0; x < nxf; x++)
        {
            for(integer z=nz0; z < nzf; z++)
            {
                const real lrho = rho_BL(rho_ref, z, x, y, dimmz, dimmx);
                
                const real stx  = stencil_X( SX, sxptr, dxi, z, x, y, dimmz, dimmx);
                const real sty  = stencil_Y( SY, syptr, dyi, z, x, y, dimmz, dimmx);
                const real stz  = stencil_Z( SZ, szptr, dzi, z, x, y, dimmz, dimmx);
                
                v_ref.bl.u[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
    ///////////////////////////////////////

#if defined(_OPENACC)
    const integer start  = ((nzf-nz0) + 2*HALO) * ((nxf-nx0) + 2*HALO) * (ny0 - HALO);
    const integer end    = ((nzf-nz0) + 2*HALO) * ((nxf-nx0) + 2*HALO) * (nyf + HALO);
    const integer nelems = end - start;

    #pragma acc data copy(calculated[start:nelems]) \
                     copyin(szptr[start:nelems], sxptr[start:nelems], syptr[start:nelems], rho[start:nelems])
    {
#endif
        compute_component_vcell_BL( v_cal.bl.u, szptr, sxptr, syptr, rho_ref, 
            dt, dzi, dxi, dyi, 
            nz0, nzf, nx0, nxf, ny0, nyf, 
            SZ, SX, SY, dimmz, dimmx, phase);
#if defined(_OPENACC)
    }
#endif

    TEST_ASSERT_EQUAL_FLOAT_ARRAY( v_ref.bl.u, v_cal.bl.u, nelems );
}

TEST(propagator, velocity_propagator_ignored_on_purpose)
{
    const real     dt  = 1.0;
    const real     dzi = 1.0;
    const real     dxi = 1.0;
    const real     dyi = 1.0;
    const integer  nz0 = HALO;
    const integer  nzf = dimmz-HALO;
    const integer  nx0 = HALO;
    const integer  nxf = dimmx-HALO;
    const integer  ny0 = HALO;
    const integer  nyf = dimmy-HALO;
    const phase_t  phase = TWO;

    //TODO: implement reference part

    velocity_propagator(v_cal, s_ref, c_ref, rho_ref, 
            dt, dzi, dxi, dyi, 
            nz0, nzf, nx0, nxf, ny0, nyf, 
            dimmz, dimmx, phase);

    TEST_IGNORE();
}

TEST(propagator, stress_update)
{   
    const real dt = 1.0;
    const real c1 = 1.0;
    const real c2 = 2.0;
    const real c3 = 3.0;
    const real c4 = 4.0;
    const real c5 = 6.0;
    const real c6 = 6.0;

    const real u_x = 5.0;
    const real v_x = 6.0;
    const real w_x = 7.0;
    
    const real u_y = 8.0;
    const real v_y = 9.0;
    const real w_y = 10.0;
    
    const real u_z = 11.0;
    const real v_z = 12.0;
    const real w_z = 13.0;

    // REFERENCE CALCULATION -DON'T TOUCH-
    for (integer y = 0; y < dimmy-1; y++)
    for (integer x = 0; x < dimmx-1; x++)
    for (integer z = 0; z < dimmz-1; z++)
    {
        const real accum = (dt * c1 * u_x) +
                           (dt * c2 * v_y) +
                           (dt * c3 * w_z) +
                           (dt * c4 * (w_y + v_z)) +
                           (dt * c5 * (w_x + u_z)) +
                           (dt * c6 * (v_x + u_y));
        s_ref.tr.xz[IDX(z,x,y,dimmz,dimmx)] += accum;
    }
    ////////////////////////////////////
   
    for (integer y = 0; y < dimmy-1; y++)
    for (integer x = 0; x < dimmx-1; x++)
    for (integer z = 0; z < dimmz-1; z++)
    {
        stress_update( s_cal.tr.xz, c1, c2, c3, c4, c5, c6, 
                z, x, y, dt, u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z, 
                dimmz, dimmx );
    }

    TEST_ASSERT_EQUAL_FLOAT_ARRAY( s_ref.tr.xz, s_cal.tr.xz, nelems );
}

TEST(propagator, stress_propagator_ignored_on_purpose)
{
    TEST_IGNORE_MESSAGE("Implement this test in the future.");
}

TEST(propagator, cell_coeff_BR)
{
    for (integer y = 0; y < dimmy-1; y++)
    for (integer x = 0; x < dimmx-1; x++)
    for (integer z = 0; z < dimmz-1; z++)
    {
        const real ref = ( 1.0f / (2.5f *( c_ref.c11[IDX(z  ,x  ,y  ,dimmz,dimmx)] +
                                           c_ref.c11[IDX(z  ,x+1,y  ,dimmz,dimmx)] +
                                           c_ref.c11[IDX(z+1,x  ,y  ,dimmz,dimmx)] +
                                           c_ref.c11[IDX(z+1,x+1,y  ,dimmz,dimmx)])) );
        
        const real cal = cell_coeff_BR( c_ref.c11, z, x, y, dimmz, dimmx);

        TEST_ASSERT_EQUAL_FLOAT(ref, cal);
    }
}

TEST(propagator, cell_coeff_TL)
{
    for (integer y = 0; y < dimmy-1; y++)
    for (integer x = 0; x < dimmx-1; x++)
    for (integer z = 0; z < dimmz-1; z++)
    {
        const real ref = ( 1.0f / c_ref.c11[IDX(z,x,y,dimmz,dimmx)]);

        const real cal = cell_coeff_TL( c_ref.c11, z, x, y, dimmz, dimmx);

        TEST_ASSERT_EQUAL_FLOAT(ref, cal);
    }
}

TEST(propagator, cell_coeff_BL)
{
    for (integer y = 0; y < dimmy-1; y++)
    for (integer x = 0; x < dimmx-1; x++)
    for (integer z = 0; z < dimmz-1; z++)
    {
        const real ref = ( 1.0f / (2.5f *( c_ref.c11[IDX(z  ,x  ,y  ,dimmz,dimmx)] +
                                           c_ref.c11[IDX(z  ,x  ,y+1,dimmz,dimmx)] +
                                           c_ref.c11[IDX(z+1,x  ,y  ,dimmz,dimmx)] +
                                           c_ref.c11[IDX(z+1,x  ,y+1,dimmz,dimmx)])) );

        const real cal = cell_coeff_BL( c_ref.c11, z, x, y, dimmz, dimmx);

        TEST_ASSERT_EQUAL_FLOAT(ref, cal);
    }
}

TEST(propagator, cell_coeff_TR)
{
    for (integer y = 0; y < dimmy-1; y++)
    for (integer x = 0; x < dimmx-1; x++)
    for (integer z = 0; z < dimmz-1; z++)
    {
        const real ref = ( 1.0f / (2.5f *( c_ref.c11[IDX(z  ,x  ,y  ,dimmz,dimmx)] +
                                           c_ref.c11[IDX(z  ,x+1,y  ,dimmz,dimmx)] +
                                           c_ref.c11[IDX(z  ,x  ,y+1,dimmz,dimmx)] +
                                           c_ref.c11[IDX(z  ,x+1,y+1,dimmz,dimmx)])) );

        const real cal = cell_coeff_TR( c_ref.c11, z, x, y, dimmz, dimmx);

        TEST_ASSERT_EQUAL_FLOAT(ref, cal);
    }
}

TEST(propagator, cell_coeff_ARTM_BR)
{
    for (integer y = 0; y < dimmy-1; y++)
    for (integer x = 0; x < dimmx-1; x++)
    for (integer z = 0; z < dimmz-1; z++)
    {
        const real ref = ((1.0f / c_ref.c11[IDX(z  ,x  ,y  ,dimmz,dimmx)] +
                           1.0f / c_ref.c11[IDX(z  ,x+1,y  ,dimmz,dimmx)] +
                           1.0f / c_ref.c11[IDX(z+1,x  ,y  ,dimmz,dimmx)] +
                           1.0f / c_ref.c11[IDX(z+1,x+1,y  ,dimmz,dimmx)])* 0.25f);

        const real cal = cell_coeff_ARTM_BR( c_ref.c11, z, x, y, dimmz, dimmx);

        TEST_ASSERT_EQUAL_FLOAT(ref, cal);
    }
}

TEST(propagator, cell_coeff_ARTM_TL)
{
    for (integer y = 0; y < dimmy-1; y++)
    for (integer x = 0; x < dimmx-1; x++)
    for (integer z = 0; z < dimmz-1; z++)
    {
        const real ref = (1.0f / c_ref.c11[IDX(z,x,y,dimmz,dimmx)]);

        const real cal = cell_coeff_ARTM_TL( c_ref.c11, z, x, y, dimmz, dimmx);

        TEST_ASSERT_EQUAL_FLOAT(ref, cal);
    }
}

TEST(propagator, cell_coeff_ARTM_BL)
{
    for (integer y = 0; y < dimmy-1; y++)
    for (integer x = 0; x < dimmx-1; x++)
    for (integer z = 0; z < dimmz-1; z++)
    {
        const real ref = ((1.0f / c_ref.c11[IDX(z  ,x  ,y  ,dimmz,dimmx)] +
                           1.0f / c_ref.c11[IDX(z  ,x  ,y+1,dimmz,dimmx)] +
                           1.0f / c_ref.c11[IDX(z+1,x  ,y  ,dimmz,dimmx)] +
                           1.0f / c_ref.c11[IDX(z+1,x  ,y+1,dimmz,dimmx)])* 0.25f);

        const real cal = cell_coeff_ARTM_BL( c_ref.c11, z, x, y, dimmz, dimmx);

        TEST_ASSERT_EQUAL_FLOAT(ref, cal);
    }
}

TEST(propagator, cell_coeff_ARTM_TR)
{
    for (integer y = 0; y < dimmy-1; y++)
    for (integer x = 0; x < dimmx-1; x++)
    for (integer z = 0; z < dimmz-1; z++)
    {
        const real ref = ((1.0f / c_ref.c11[IDX(z  ,x  ,y  ,dimmz,dimmx)] +
                           1.0f / c_ref.c11[IDX(z  ,x+1,y  ,dimmz,dimmx)] +
                           1.0f / c_ref.c11[IDX(z  ,x  ,y+1,dimmz,dimmx)] +
                           1.0f / c_ref.c11[IDX(z  ,x+1,y+1,dimmz,dimmx)])* 0.25f);

        const real cal = cell_coeff_ARTM_TR( c_ref.c11, z, x, y, dimmz, dimmx);

        TEST_ASSERT_EQUAL_FLOAT(ref, cal);
    }
}


////// TESTS RUNNER //////
TEST_GROUP_RUNNER(propagator)
{
    RUN_TEST_CASE(propagator, IDX);

    RUN_TEST_CASE(propagator, stencil_Z);
    RUN_TEST_CASE(propagator, stencil_X);
    RUN_TEST_CASE(propagator, stencil_Y);

    RUN_TEST_CASE(propagator, rho_BL);
    RUN_TEST_CASE(propagator, rho_TR);
    RUN_TEST_CASE(propagator, rho_BR);
    RUN_TEST_CASE(propagator, rho_TL);

    RUN_TEST_CASE(propagator, compute_component_vcell_TL);
    RUN_TEST_CASE(propagator, compute_component_vcell_TR);
    RUN_TEST_CASE(propagator, compute_component_vcell_BR);
    RUN_TEST_CASE(propagator, compute_component_vcell_BL);

    RUN_TEST_CASE(propagator, velocity_propagator_ignored_on_purpose);
    
    RUN_TEST_CASE(propagator, stress_update);

    RUN_TEST_CASE(propagator, stress_propagator_ignored_on_purpose);

    RUN_TEST_CASE(propagator, cell_coeff_BR);
    RUN_TEST_CASE(propagator, cell_coeff_TL);
    RUN_TEST_CASE(propagator, cell_coeff_BL);
    RUN_TEST_CASE(propagator, cell_coeff_TR);

    RUN_TEST_CASE(propagator, cell_coeff_ARTM_BR);
    RUN_TEST_CASE(propagator, cell_coeff_ARTM_TL);
    RUN_TEST_CASE(propagator, cell_coeff_ARTM_BL);
    RUN_TEST_CASE(propagator, cell_coeff_ARTM_TR);
}
