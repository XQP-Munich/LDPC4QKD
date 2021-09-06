//
// Created by alice on, 07.05.21.
//

// Google Test framework
#include <gtest/gtest.h>
#include "helpers_for_testing.hpp"

// Standard library
#include <iostream>

// To be tested
#include "read_scsmat_format.hpp"
#include "fortest_autogen_ldpc_matrix_csc.hpp"

using namespace HelpersForTests;
using namespace LDPC4QKD;


TEST(test_read_scmat_format, read_matrix_from_cscmat) {
    auto pair = read_matrix_from_cscmat("./LDPC_code_for_testing_2048x6144.cscmat");
    auto colptr = pair.first;
    auto row_idx = pair.second;

    EXPECT_EQ(std::vector<decltype(colptr)::value_type>(AutogenLDPC::colptr.begin(), AutogenLDPC::colptr.end()),
            colptr);

    EXPECT_EQ(std::vector<decltype(row_idx)::value_type>(AutogenLDPC::row_idx.begin(), AutogenLDPC::row_idx.end()),
              row_idx);
}
