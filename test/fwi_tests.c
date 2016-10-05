#include "fwi_tests.h"


/// HELPER FUNCTIONS /////
bool    is_negative(Float_t f)  { return f.i < 0;               }
int32_t raw_mantissa(Float_t f) { return f.i & ((1 << 23) - 1); }
int32_t raw_exponent(Float_t f) { return (f.i >> 23) & 0xFF;    }
         
bool assert_float_equal_ULPs(float A, float B, int maxUlpsDiff)
{
    Float_t uA;
    Float_t uB;
    
    uA.f = A;
    uB.f = B;

    // NaN or Infinite means they do not match.
    if (isnan(A) || isnan(B) || (A == INFINITY) || (B == INFINITY))
        return false;

    // Different signs means they do not match.
    if (is_negative(uA) != is_negative(uB))
    {
        // Check for equality to make sure +0==-0
        if (A == B)
            return true;
        return false;
    }
                                                      
    // Find the difference in ULPs.
    int ulpsDiff = abs(uA.i - uB.i);
    if (ulpsDiff <= maxUlpsDiff)
        return true;
                                                       
    return false;
}

int diff_in_ULPs(float A, float B)
{
    Float_t uA, uB;
    uA.f = A;
    uB.f = B;

    return abs(uA.i - uB.i);
}

void assert_equal_float_array( float* ref, float* opt, const int nelems, const char* file, const int line )
{
    // OBS: the maximum ULP difference observed is 16384 and it seems to come from difference
    //      optimizations applied by NVCC compiler in "stress_update" function in scell_BR
    //      causes an accumulation of errors.
    int maxULP = 16384;
    int max_diff = 0;

    for (int elem = 0; elem < nelems; elem++)
    {
        float ref_e = ref[elem];
        float opt_e = opt[elem];
        const int diff = diff_in_ULPs(ref_e, opt_e);

        if (max_diff < diff) max_diff = diff;

        if ( !assert_float_equal_ULPs(ref_e, opt_e, maxULP) )
        {
            fprintf(stderr, "ERROR:%s:%d: in element %d : ref %e !=\t %e \topt by %d ULP\n",
                   file, line, elem, ref_e, opt_e, diff );

            Float_t uref, uopt;
            uref.f = ref_e;
            uopt.f = opt_e;
            fprintf(stderr, "REF: sign %d mantissa %023d exponent %08d\n",
                is_negative(uref), raw_mantissa(uref), raw_exponent(uref));
            fprintf(stderr, "OPT: sign %d mantissa %023d exponent %08d\n",
                is_negative(uopt), raw_mantissa(uopt), raw_exponent(uopt));


            exit(1); //UNITY_FAIL_AND_BAIL;
        }
    }
    fprintf(stdout, "\nMAX DIFF: %d ULPs", max_diff);
}

///// UNITY TEST RUNNER //////
static void run_all_tests(void)
{
    RUN_TEST_GROUP(common);
    RUN_TEST_GROUP(propagator);
    RUN_TEST_GROUP(kernel);
}

int main(int argc, const char* argv[])
{
    return UnityMain(argc, argv, run_all_tests);
}
