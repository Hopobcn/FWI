#ifndef _FWI_TESTS_H_
#define _FWI_TESTS_H_

#include "fwi_propagator.h"
#include "unity_config.h"

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


#endif /* end of _FWI_TESTS_H_ definiton */
