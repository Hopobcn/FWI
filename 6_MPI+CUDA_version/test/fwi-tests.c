#include <unity.h>

#include "fwi_common.h"


void test_dummy_PASS()
{
    // DUMMY test
    TEST_ASSERT_EQUAL(1, 1);
}

void test_dummy_FAILS()
{
    // DUMMY test that always fails (Control purposes)
    TEST_ASSERT_EQUAL(0, 1);
}

int main()
{
    UNITY_BEGIN();

    RUN_TEST(test_dummy_PASS);
    RUN_TEST(test_dummy_FAILS);
    
    return UNITY_END();
}
