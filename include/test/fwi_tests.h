#ifndef _FWI_TESTS_H_
#define _FWI_TESTS_H_

#include "fwi_propagator.h"
#include "unity_config.h"

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
