// Google Test framework
#include <gtest/gtest.h>
#include "helpers_for_testing.hpp"

#include <iostream>

using namespace HelpersForTests;


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}


TEST(TestUtils, valgrind_check) {
    // Running this test is only useful with valgrind memory check tool enabled.
    // This test mishandles memory on purpose to make sure if the error is detected by valgrind.
    int x, y = 0;
    if (x == 42) y = 2; // Should give UninitCondition warning!
    else y++;
}


TEST(TestUtils, hash_vector) {
    std::vector<std::uint32_t> test_vec{0, 1, 2, 3, 4};
    EXPECT_EQ(hash_vector(test_vec), 3632105860);
}

TEST(TestUtils, gen_bitstring) {
    std::vector<std::uint32_t> test_vec = get_bitstring<std::uint32_t>(1234);
    EXPECT_EQ(hash_vector(test_vec), 3900352086);

    std::vector<bool> test_vec2 = get_bitstring<bool>(1234);
    EXPECT_EQ(hash_vector(test_vec2), 3900352086);
}
