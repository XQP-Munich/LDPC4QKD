//
// Created by alice on, 07.05.21.
//

// Google Test framework
#include <gtest/gtest.h>
#include "helpers_for_testing.hpp"

// Standard library
#include <iostream>

// To be tested
#include "rate_adaptive_code.hpp"

// Test cases test against constants known to be correct for the LDPC-matrix defined here:
#include "fortest_autogen_ldpc_matrix_csc.hpp"

using namespace HelpersForTests;
using namespace LDPC4QKD;

namespace {

    RateAdaptiveCode<Bit> get_code_big() {
        std::vector<std::uint32_t> colptr(AutogenLDPC::colptr.begin(), AutogenLDPC::colptr.end());
        std::vector<std::uint16_t> row_idx(AutogenLDPC::row_idx.begin(), AutogenLDPC::row_idx.end());
        return RateAdaptiveCode<Bit>(colptr, row_idx);
    }


    RateAdaptiveCode<Bit> get_code_small() {
        //    H =  [1 0 1 0 1 0 1
        //			0 1 1 0 0 1 1
        //			0 0 0 1 1 1 1]
        std::vector<std::uint32_t> colptr{0, 1, 2, 4, 5, 7, 9, 12};
        std::vector<std::uint16_t> row_idx{0, 1, 0, 1, 2, 0, 2, 1, 2, 0, 1, 2};
        return RateAdaptiveCode<Bit>(colptr, row_idx);
    }

}

//TEST(rate_adaptive_code, TMPTMPTMPTMTPTMP) { // this test accesses private fields.
//    auto H = get_code_big();
//    EXPECT_EQ(hash_vector(H.colptr), 736283749);
//    EXPECT_EQ(hash_vector(H.row_idx), 4281948431);
//}


TEST(rate_adaptive_code, decode_test_small) {
    auto H = get_code_small();

    std::vector<Bit> x{1, 1, 1, 1, 0, 0, 0}; // true data to be sent
    std::vector<Bit> syndrome;
    H.encode_no_ra(x, syndrome);

    std::vector<Bit> x_noised{1, 1, 1, 1, 0, 0, 1}; // distorted data
    double p = 1. / 7; // channel error probability (we flipped 1 symbol out of 7)

    double vlog = log((1 - p) / p);
    std::vector<double> llrs(x.size());
    for (std::size_t i{}; i < llrs.size(); ++i) {
        llrs[i] = vlog * (1 - 2 * x_noised[i]); // log likelihood ratios
    }

    std::vector<Bit> solution;
    bool success = H.decode_at_current_rate(llrs, syndrome, solution);
    EXPECT_TRUE(success);
    EXPECT_EQ(solution, x);
}

TEST(rate_adaptive_code, decode_test_big) {
    auto H = get_code_big();

    std::vector<Bit> x = get_bitstring(H.getNCols()); // true data to be sent
    std::vector<Bit> syndrome;
    H.encode_no_ra(x, syndrome);

    constexpr double p = 0.04; // channel error probability
    std::vector<Bit> x_noised = x; // distorted data
    noise_bitstring_inplace(x_noised, p);

    double vlog = log((1 - p) / p);
    std::vector<double> llrs(x.size());
    for (std::size_t i{}; i < llrs.size(); ++i) {
        llrs[i] = vlog * (1 - 2 * x_noised[i]); // log likelihood ratios
    }

    std::vector<Bit> solution;
    bool success = H.decode_at_current_rate(llrs, syndrome, solution);
    EXPECT_TRUE(success);
    EXPECT_EQ(solution, x);
    for (std::size_t i{}; i < x.size(); ++i) {
        if (solution[i] != x[i]) {
            std::cout << "Error at bit position " << i << std::endl;
        }
    }
}


TEST(rate_adaptive_code, encode_no_ra) {
    auto H = get_code_big();
    std::vector<Bit> in = get_bitstring(H.getNCols());
    std::vector<Bit> out(H.get_NRows_mother_matrix());

    std::cout << hash_vector(in) << std::endl;
    H.encode_no_ra(in, out);

    EXPECT_EQ(hash_vector(out), 2814594723);
}


TEST(rate_adaptive_code, encode_current_rate) {
    auto H = get_code_big();
    std::vector<Bit> in = get_bitstring(H.getNCols());
    std::vector<Bit> out(H.get_NRows_mother_matrix());

    H.encode_at_current_rate(in, out);

    EXPECT_EQ(hash_vector(out), 2814594723);

    // TODO change rate!!
}


TEST(rate_adaptive_code, init_pos_CN_pos_VN) {
    auto H = get_code_small();

    std::vector<std::vector<decltype(H)::MatrixIndex>> expect_posCN{{0},
                                                                    {1},
                                                                    {0, 1},
                                                                    {2},
                                                                    {0, 2},
                                                                    {1, 2},
                                                                    {0, 1, 2}};
    std::vector<std::vector<decltype(H)::MatrixIndex>> expect_posVN{{0, 2, 4, 6},
                                                                    {1, 2, 5, 6},
                                                                    {3, 4, 5, 6}};
    EXPECT_EQ(H.getPosVarn(), expect_posVN);
    EXPECT_EQ(H.getPosCheckn(), expect_posCN);
}


TEST(rate_adaptive_code, getters) {
    auto H = get_code_big();
    EXPECT_EQ(H.get_NRows_mother_matrix(), 2048);
    EXPECT_EQ(H.getNCols(), 6144);

//    H.get_current_n_rate_adapted_rows();
}

TEST(rate_adaptive_code, encode_with_ra) {
    auto H = get_code_big();

    std::vector<Bit> input = get_bitstring(H.getNCols()); // true data to be sent

    // storage for syndrome. Initialize with arbitrary values, which must be overwritten by encoder.
    std::vector<Bit> syndrome = get_bitstring(H.get_NRows_mother_matrix());

    // test agreement with encoder that doesn't use rate adaption
    H.encode_with_ra(input, syndrome, H.get_NRows_mother_matrix());

    EXPECT_EQ(hash_vector(syndrome), 2814594723);

    // check that invalid requests lead to exceptions
    EXPECT_ANY_THROW(H.encode_with_ra({true, false}, syndrome, -1));  // invalid input size
    EXPECT_ANY_THROW(H.encode_with_ra(input, syndrome, -1)); // invalid requested size
    EXPECT_ANY_THROW(H.encode_with_ra(input, syndrome, H.get_NRows_mother_matrix() + 1)); // too big requested size

    H.encode_with_ra(input, syndrome, H.get_NRows_mother_matrix() / 2 + 1);
    EXPECT_EQ(hash_vector(syndrome), 2814594723);

}


TEST(rate_adaptive_code, rate_adaption) {
    // TODO
    FAIL();
//    auto H = get_code_big(give rate adaption initial);
    // check sizes
    // change rate adaption
    // check sizes

}
