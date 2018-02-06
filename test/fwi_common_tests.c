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
