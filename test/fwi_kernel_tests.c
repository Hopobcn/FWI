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
#include "test/fwi_tests.h"

#include <unity.h>
#include <unity_fixture.h>

#include "fwi/fwi_kernel.h"



TEST_GROUP(kernel);

TEST_SETUP(kernel)
{
    nelems = dimmz * dimmx * dimmy;

    alloc_memory_shot(dimmz, dimmx, dimmy, &c_ref, &s_ref, &v_ref, &rho_ref);
    alloc_memory_shot(dimmz, dimmx, dimmy, &c_cal, &s_cal, &v_cal, &rho_cal);

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

TEST_TEAR_DOWN(kernel)
{
    free_memory_shot(&c_ref, &s_ref, &v_ref, &rho_ref);
    free_memory_shot(&c_cal, &s_cal, &v_cal, &rho_cal);
}

TEST(kernel, set_array_to_random_real)
{
    #define NELEMS 100
    real array_ref[NELEMS];
    real array_cal[NELEMS];

    srand(0);
    {
        const real randvalue = rand() / (1.0 * RAND_MAX);

        for( integer i = 0; i < NELEMS; i++ )
            array_ref[i] = randvalue;
    }

    srand(0);
#if defined(_OPENACC)
    #pragma acc data copy(array_cal[:NELEMS])
#endif
    {
        set_array_to_random_real( array_cal, NELEMS );
    }

    TEST_ASSERT_EQUAL_FLOAT_ARRAY( array_ref, array_cal, NELEMS );
}

TEST(kernel, set_array_to_constant)
{
    #define NELEMS 100
    real array_ref[NELEMS];
    real array_cal[NELEMS];

    const real value = 3.0;

    {
        for( integer i = 0; i < NELEMS; i++ )
            array_ref[i] = value;
    }

#if defined(_OPENACC)
    #pragma acc data copy(array_cal[:NELEMS])
#endif
    {
        set_array_to_constant( array_cal,value, NELEMS );
    }

    TEST_ASSERT_EQUAL_FLOAT_ARRAY( array_ref, array_cal, NELEMS );
}

////// TESTS RUNNER //////
TEST_GROUP_RUNNER(kernel)
{
    RUN_TEST_CASE(kernel, set_array_to_random_real);
    RUN_TEST_CASE(kernel, set_array_to_constant);
}
