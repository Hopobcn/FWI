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

#ifndef _FWI_TESTS_H_
#define _FWI_TESTS_H_

#include "fwi/fwi_propagator.h"
#include "test/unity_config.h"

#include <unity.h>
#include <unity_fixture.h>

#include <stdbool.h>


// global vars
extern const integer dimmz;
extern const integer dimmx;
extern const integer dimmy;
extern       integer nelems;

extern v_t v_ref;
extern s_t s_ref;
extern coeff_t c_ref;
extern real* rho_ref;

extern v_t v_cal;
extern s_t s_cal;
extern coeff_t c_cal;
extern real* rho_cal;

// Helper functions
void init_array( real* restrict array, const integer length );
void copy_array( real* restrict dest, real* restrict src, const integer length );


typedef union
{
    int32_t i;
    float   f;
#ifdef DEBUG
    struct
    {  // Bitfields for exploration.
        uint32_t mantissa : 23;
        uint32_t exponent : 8;
        uint32_t sign : 1;
    } parts;
#endif
} Float_t;

#define CUSTOM_ASSERT_EQUAL_FLOAT_ARRAY(ref,opt,nelems)                     \
{                                                                           \
    assert_equal_float_array( (ref), (opt), (nelems), __FILE__, __LINE__ ); \
}

void assert_equal_float_array( float* ref, float* opt, const int nelems, const char* file, const int line );


#endif /* end of _FWI_TESTS_H_ definiton */
