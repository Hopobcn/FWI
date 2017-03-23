#include <unity.h>
#include <unity_fixture.h>
#include <stdbool.h>

#include "fwi/fwi_common.h"

TEST_GROUP(common);


TEST_SETUP(common)
{
}

TEST_TEAR_DOWN(common)
{
}

TEST(common, roundup)
{
    TEST_ASSERT_EQUAL_INT(0,  roundup(0, 0));
  //TEST_ASSERT_EQUAL_INT(32, roundup(0, 32));
    TEST_ASSERT_EQUAL_INT(32, roundup(31, 32));
    TEST_ASSERT_EQUAL_INT(32, roundup(32, 32));
    TEST_ASSERT_EQUAL_INT(64, roundup(33, 32));

    TEST_ASSERT_EQUAL_INT(  HALO, roundup(HALO-1, HALO));
    TEST_ASSERT_EQUAL_INT(  HALO, roundup(HALO,   HALO));
    TEST_ASSERT_EQUAL_INT(2*HALO, roundup(HALO+1, HALO));
}

bool is_aligned(void* p, size_t n)
{
    return (size_t)p % n == 0;
}

TEST(common, aligned_allocations)
{
    extent_t req = make_extent(284,284,284);
    dim_t dim;
    const size_t aligment = 128;

    float* h_ptr = (float*) malloc3d_host(&dim, aligment, HALO, req);

    float* row_0 = &h_ptr[IDX(HALO,0,0,dim)];
    float* row_1 = &h_ptr[IDX(HALO,1,0,dim)];
    float* row_2 = &h_ptr[IDX(HALO,2,0,dim)];
    float* row_3 = &h_ptr[IDX(HALO,3,0,dim)];

    TEST_ASSERT_TRUE( is_aligned(row_0, aligment) );
    TEST_ASSERT_TRUE( is_aligned(row_1, aligment) );
    TEST_ASSERT_TRUE( is_aligned(row_2, aligment) );
    TEST_ASSERT_TRUE( is_aligned(row_3, aligment) );


    float* d_ptr = (float*) malloc3d_device(&dim, aligment, HALO, req, (void*)h_ptr);

    float* row_4 = &d_ptr[IDX(HALO,0,0,dim)];
    float* row_5 = &d_ptr[IDX(HALO,1,0,dim)];
    float* row_6 = &d_ptr[IDX(HALO,2,0,dim)];
    float* row_7 = &d_ptr[IDX(HALO,3,0,dim)];

    TEST_ASSERT_TRUE( is_aligned(row_4, aligment) );
    TEST_ASSERT_TRUE( is_aligned(row_5, aligment) );
    TEST_ASSERT_TRUE( is_aligned(row_6, aligment) );
    TEST_ASSERT_TRUE( is_aligned(row_7, aligment) );

    free3d_device(h_ptr);

    free3d_host(h_ptr);
}

TEST_GROUP_RUNNER(common)
{
    RUN_TEST_CASE(common, roundup);
    RUN_TEST_CASE(common, aligned_allocations);
}
