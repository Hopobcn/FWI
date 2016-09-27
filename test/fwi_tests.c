#include <unity.h>
#include <unity_fixture.h>

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
