#include <unity.h>

#include "fwi_common.h"

// Function Declarations (Prototypes)
void test_dummy_PASS(void);
void test_dummy_FAILS(void);


// Function Definitions 
void test_dummy_PASS(void)
{
    // DUMMY test
    TEST_ASSERT_EQUAL(1, 1);
}

void test_dummy_FAILS(void)
{
    // DUMMY test that always fails (Control purposes)
    TEST_ASSERT_EQUAL(0, 1);
}

int main(int argc, char* argv[])
{
    UNITY_BEGIN();

    RUN_TEST(test_dummy_PASS);
    RUN_TEST(test_dummy_FAILS);
    
    return UNITY_END();
}
