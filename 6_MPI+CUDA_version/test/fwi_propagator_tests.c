#include <unity.h>
#include <unity_fixture.h>

#include "fwi_propagator.h"

TEST_GROUP(propagator);


////// SetUp/SetDown /////
const integer dimmz = 32;
const integer dimmx = 16;
const integer dimmy = 16;
      integer nelems;

real* calculated;
real* vptr;
real* szptr;
real* sxptr;
real* syptr;
real* rho;

TEST_SETUP(propagator)
{
    nelems = dimmz * dimmx * dimmy;

    calculated = (real*) malloc( nelems * sizeof(real) );
    vptr  = (real*) malloc( nelems * sizeof(real) );
    szptr = (real*) malloc( nelems * sizeof(real) ); 
    sxptr = (real*) malloc( nelems * sizeof(real) ); 
    syptr = (real*) malloc( nelems * sizeof(real) ); 
    rho   = (real*) malloc( nelems * sizeof(real) );

    for (int i = 0; i < nelems; i++) calculated[i] = i/(1.0+i);
    for (int i = 0; i < nelems; i++) vptr[i] = calculated[i];
    for (int i = 0; i < nelems; i++) szptr[i] = calculated[i];
    for (int i = 0; i < nelems; i++) sxptr[i] = calculated[i];
    for (int i = 0; i < nelems; i++) syptr[i] = calculated[i];
    for (int i = 0; i < nelems; i++) rho[i] = calculated[i];
}

TEST_TEAR_DOWN(propagator)
{
    free(calculated);
    free(vptr);
    free(szptr);
    free(sxptr);
    free(syptr);
    free(rho);
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
        TEST_ASSERT_EQUAL_FLOAT(stencil_Z(FORWARD,  calculated, dzi, z,   x, y, dimmz, dimmx),
                                stencil_Z(BACKWARD, calculated, dzi, z+1, x, y, dimmz, dimmx));
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
        TEST_ASSERT_EQUAL_FLOAT(stencil_X(FORWARD,  calculated, dzi, z, x,   y, dimmz, dimmx),
                                stencil_X(BACKWARD, calculated, dzi, z, x+1, y, dimmz, dimmx));
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
        TEST_ASSERT_EQUAL_FLOAT(stencil_Y(FORWARD,  calculated, dzi, z, x, y,   dimmz, dimmx),
                                stencil_Y(BACKWARD, calculated, dzi, z, x, y+1, dimmz, dimmx));
    }
}

TEST(propagator, rho_BL)
{
    for (integer y = 0; y < dimmy;   y++)
    for (integer x = 0; x < dimmx;   x++)
    for (integer z = 0; z < dimmz-1; z++)
    {
        TEST_ASSERT_EQUAL_FLOAT((2.0f/(calculated[IDX(z,x,y,dimmz,dimmx)]+calculated[IDX(z+1,x,y,dimmz,dimmx)])),
                                rho_BL(calculated, z, x, y, dimmz, dimmx));
    }
}

TEST(propagator, rho_TR)
{
    for (integer y = 0; y < dimmy;   y++)
    for (integer x = 0; x < dimmx-1; x++)
    for (integer z = 0; z < dimmz;   z++)
    {
        TEST_ASSERT_EQUAL_FLOAT((2.0f/(calculated[IDX(z,x,y,dimmz,dimmx)]+calculated[IDX(z,x+1,y,dimmz,dimmx)])),
                                rho_TR(calculated, z, x, y, dimmz, dimmx));
    }
}

TEST(propagator, rho_BR)
{
    for (integer y = 0; y < dimmy-1; y++)
    for (integer x = 0; x < dimmx-1; x++)
    for (integer z = 0; z < dimmz-1; z++)
    {
        const real ref = ( 8.0f / ( calculated[IDX(z  ,x  ,y  ,dimmz,dimmx)] +
                                    calculated[IDX(z+1,x  ,y  ,dimmz,dimmx)] +
                                    calculated[IDX(z  ,x+1,y  ,dimmz,dimmx)] +
                                    calculated[IDX(z  ,x  ,y+1,dimmz,dimmx)] +
                                    calculated[IDX(z  ,x+1,y+1,dimmz,dimmx)] +
                                    calculated[IDX(z+1,x+1,y  ,dimmz,dimmx)] +
                                    calculated[IDX(z+1,x  ,y+1,dimmz,dimmx)] +
                                    calculated[IDX(z+1,x+1,y+1,dimmz,dimmx)]) );
        
        const real cal = rho_BR(calculated, z, x, y, dimmz, dimmx);

        TEST_ASSERT_EQUAL_FLOAT(ref, cal);
    }
}

TEST(propagator, rho_TL)
{
    for (integer y = 0; y < dimmy-1; y++)
    for (integer x = 0; x < dimmx;   x++)
    for (integer z = 0; z < dimmz;   z++)
    {
        TEST_ASSERT_EQUAL_FLOAT((2.0f/(calculated[IDX(z,x,y,dimmz,dimmx)]+calculated[IDX(z,x,y+1,dimmz,dimmx)])),
                                rho_TL(calculated, z, x, y, dimmz, dimmx));
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

    // REFERENCE CALCULATION -DON'T TOUCH-
    for(integer y=ny0; y < nyf; y++)
    {
        for(integer x=nx0; x < nxf; x++)
        {
            for(integer z=nz0; z < nzf; z++)
            {
                const real lrho = rho_TL(rho, z, x, y, dimmz, dimmx);
                
                const real stx  = stencil_X( SX, sxptr, dxi, z, x, y, dimmz, dimmx);
                const real sty  = stencil_Y( SY, syptr, dyi, z, x, y, dimmz, dimmx);
                const real stz  = stencil_Z( SZ, szptr, dzi, z, x, y, dimmz, dimmx);
                
                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
    ///////////////////////////////////////

    compute_component_vcell_TL( calculated, szptr, sxptr, syptr, rho, 
        dt, dzi, dxi, dyi, 
        nz0, nzf, nx0, nxf, ny0, nyf, 
        SZ, SX, SY, dimmz, dimmx, phase);

    TEST_ASSERT_EQUAL_FLOAT_ARRAY( vptr, calculated, nelems );
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

    // REFERENCE CALCULATION -DON'T TOUCH-
    for(integer y=ny0; y < nyf; y++)
    {
        for(integer x=nx0; x < nxf; x++)
        {
            for(integer z=nz0; z < nzf; z++)
            {
                const real lrho = rho_TR(rho, z, x, y, dimmz, dimmx);
                
                const real stx  = stencil_X( SX, sxptr, dxi, z, x, y, dimmz, dimmx);
                const real sty  = stencil_Y( SY, syptr, dyi, z, x, y, dimmz, dimmx);
                const real stz  = stencil_Z( SZ, szptr, dzi, z, x, y, dimmz, dimmx);
                
                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
    ///////////////////////////////////////

    compute_component_vcell_TR( calculated, szptr, sxptr, syptr, rho, 
        dt, dzi, dxi, dyi, 
        nz0, nzf, nx0, nxf, ny0, nyf, 
        SZ, SX, SY, dimmz, dimmx, phase);

    TEST_ASSERT_EQUAL_FLOAT_ARRAY( vptr, calculated, nelems );
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

    // REFERENCE CALCULATION -DON'T TOUCH-
    for(integer y=ny0; y < nyf; y++)
    {
        for(integer x=nx0; x < nxf; x++)
        {
            for(integer z=nz0; z < nzf; z++)
            {
                const real lrho = rho_BR(rho, z, x, y, dimmz, dimmx);
                
                const real stx  = stencil_X( SX, sxptr, dxi, z, x, y, dimmz, dimmx);
                const real sty  = stencil_Y( SY, syptr, dyi, z, x, y, dimmz, dimmx);
                const real stz  = stencil_Z( SZ, szptr, dzi, z, x, y, dimmz, dimmx);
                
                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
    ///////////////////////////////////////

    compute_component_vcell_BR( calculated, szptr, sxptr, syptr, rho, 
        dt, dzi, dxi, dyi, 
        nz0, nzf, nx0, nxf, ny0, nyf, 
        SZ, SX, SY, dimmz, dimmx, phase);

    TEST_ASSERT_EQUAL_FLOAT_ARRAY( vptr, calculated, nelems );
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

    // REFERENCE CALCULATION -DON'T TOUCH-
    for(integer y=ny0; y < nyf; y++)
    {
        for(integer x=nx0; x < nxf; x++)
        {
            for(integer z=nz0; z < nzf; z++)
            {
                const real lrho = rho_BL(rho, z, x, y, dimmz, dimmx);
                
                const real stx  = stencil_X( SX, sxptr, dxi, z, x, y, dimmz, dimmx);
                const real sty  = stencil_Y( SY, syptr, dyi, z, x, y, dimmz, dimmx);
                const real stz  = stencil_Z( SZ, szptr, dzi, z, x, y, dimmz, dimmx);
                
                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
    ///////////////////////////////////////

    compute_component_vcell_BL( calculated, szptr, sxptr, syptr, rho, 
        dt, dzi, dxi, dyi, 
        nz0, nzf, nx0, nxf, ny0, nyf, 
        SZ, SX, SY, dimmz, dimmx, phase);

    TEST_ASSERT_EQUAL_FLOAT_ARRAY( vptr, calculated, nelems );
}

TEST(propagator, velocity_propagator_ignored_on_purpose)
{
    TEST_IGNORE_MESSAGE("Implement this test in the future");
}

TEST(propagator, stress_update)
{
    for (integer y = 0; y < dimmy-1; y++)
    for (integer x = 0; x < dimmx;   x++)
    for (integer z = 0; z < dimmz;   z++)
    {
        //TEST_ASSERT_EQUAL_FLOAT((2.0f/(calculated[IDX(z,x,y,dimmz,dimmx)]+calculated[IDX(z,x,y+1,dimmz,dimmx)])),     rho_TL(calculated, z, x, y, dimmz, dimmx));
    }
    TEST_IGNORE();
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
}
