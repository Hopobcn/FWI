#include <unity.h>
#include <unity_fixture.h>

#include "fwi/fwi_common.h"

TEST_GROUP(common);


TEST_SETUP(common)
{
}

TEST_TEAR_DOWN(common)
{
}

TEST(common, checkErrors_ignored_on_purpose)
{
    TEST_IGNORE_MESSAGE("why I ignore this test..");
    //TEST_IGNORE();
}

TEST(common, safe_fopen)
{
    TEST_ASSERT_EQUAL(0, 0);
}

TEST(common, safe_fclose)
{
    TEST_ASSERT_EQUAL(0, 0);
}

TEST(common, safe_fwrite)
{
    TEST_ASSERT_EQUAL(0, 0);
}

TEST(common, safe_fread)
{
    TEST_ASSERT_EQUAL_FLOAT(3.1415, 3.1415);
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

// ...


TEST_GROUP_RUNNER(common)
{
    RUN_TEST_CASE(common, checkErrors_ignored_on_purpose);
    
    RUN_TEST_CASE(common, safe_fopen);
    RUN_TEST_CASE(common, safe_fclose);
    RUN_TEST_CASE(common, safe_fwrite);
    RUN_TEST_CASE(common, safe_fread);

    RUN_TEST_CASE(common, roundup);
}
