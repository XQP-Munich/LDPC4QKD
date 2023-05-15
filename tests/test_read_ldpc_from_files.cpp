//
// Created by alice on, 07.05.21.
//

// Google Test framework
#include <gtest/gtest.h>
#include "helpers_for_testing.hpp"

// Standard library
#include <iostream>

// To be tested
#include "read_ldpc_file_formats.hpp"
#include "fortest_autogen_ldpc_matrix_csc.hpp"
#include "benchmarks_error_rate/code_simulation_helpers.hpp"

using namespace HelpersForTests;
using namespace LDPC4QKD;

namespace {

    auto get_code_big_nora() {
        std::vector<std::uint32_t> colptr(AutogenLDPC::colptr.begin(), AutogenLDPC::colptr.end());
        std::vector<std::uint32_t> row_idx(AutogenLDPC::row_idx.begin(), AutogenLDPC::row_idx.end());
        return RateAdaptiveCode<Bit, uint32_t, uint32_t>(colptr, row_idx);
    }

}

TEST(test_read_ldpc_from_files, read_matrix_from_cscmat) {
    auto pair = read_matrix_from_cscmat<std::uint32_t, std::uint16_t>(
            "./LDPC_code_for_testing_2048x6144.cscmat");
    auto colptr = pair.first;
    auto row_idx = pair.second;

    EXPECT_EQ(std::vector<decltype(colptr)::value_type>(AutogenLDPC::colptr.begin(), AutogenLDPC::colptr.end()),
            colptr);

    EXPECT_EQ(std::vector<decltype(row_idx)::value_type>(AutogenLDPC::row_idx.begin(), AutogenLDPC::row_idx.end()),
              row_idx);
}

TEST(test_read_ldpc_from_files, read_rate_adaption_from_csv_) {
    auto rows_to_combine = read_rate_adaption_from_csv<std::size_t>(
            "./rate_adaption_2x6_block_6144_for_testing.csv");
    auto hash = hash_vector(rows_to_combine);
    EXPECT_EQ(hash, 453016743);
}


TEST(test_read_ldpc_from_files, read_bincsc_json_format) {
    auto H = LDPC4QKD::CodeSimulationHelpers::load_ldpc_from_json("./test_reading_bincscjson_format_block_6144_proto_2x6_313422410401.bincsc.json");
    auto H_old = get_code_big_nora();

    EXPECT_TRUE(H == H_old);
}
