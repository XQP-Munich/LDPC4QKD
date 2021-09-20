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
#include "fortest_autogen_rate_adaption.hpp"

using namespace HelpersForTests;
using namespace LDPC4QKD;

namespace {

    RateAdaptiveCode<Bit> get_code_big_nora() {
        std::vector<std::uint32_t> colptr(AutogenLDPC::colptr.begin(), AutogenLDPC::colptr.end());
        std::vector<std::uint16_t> row_idx(AutogenLDPC::row_idx.begin(), AutogenLDPC::row_idx.end());
        return RateAdaptiveCode<Bit>(colptr, row_idx);
    }

    RateAdaptiveCode<Bit> get_code_big_wra() {
        std::vector<std::uint32_t> colptr(AutogenLDPC::colptr.begin(), AutogenLDPC::colptr.end());
        std::vector<std::uint16_t> row_idx(AutogenLDPC::row_idx.begin(), AutogenLDPC::row_idx.end());
        std::vector<std::uint16_t> rows_to_combine(AutogenRateAdapt::rows.begin(), AutogenRateAdapt::rows.end());
        return RateAdaptiveCode<Bit>(colptr, row_idx, rows_to_combine);
    }


    RateAdaptiveCode<Bit> get_code_small() {
        //    H =  [1 0 1 0 1 0 1
        //			0 1 1 0 0 1 1
        //			0 0 0 1 1 1 1]
        std::vector<std::uint32_t> colptr{0, 1, 2, 4, 5, 7, 9, 12};
        std::vector<std::uint16_t> row_idx{0, 1, 0, 1, 2, 0, 2, 1, 2, 0, 1, 2};
        return RateAdaptiveCode<Bit>(colptr, row_idx);
    }


    double h2(double p) {
        return -p * ::log(p) - (1 - p) * log(1 - p);
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
    auto H = get_code_big_nora();

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
    auto H = get_code_big_nora();
    std::vector<Bit> in = get_bitstring(H.getNCols());
    std::vector<Bit> out(H.get_n_rows_mother_matrix());

    std::cout << hash_vector(in) << std::endl;
    H.encode_no_ra(in, out);

    EXPECT_EQ(hash_vector(out), 2814594723);
}


TEST(rate_adaptive_code, encode_current_rate) {
    auto H = get_code_big_wra();
    std::vector<Bit> in = get_bitstring(H.getNCols());
    std::vector<Bit> out(H.get_n_rows_mother_matrix());

    H.encode_at_current_rate(in, out);

    EXPECT_EQ(hash_vector(out), 2814594723);

    const auto n_line_combs = static_cast<std::size_t>(H.get_n_rows_mother_matrix() * 0.3);
    H.set_rate(n_line_combs);

    H.encode_at_current_rate(in, out);
    EXPECT_EQ(hash_vector(out), 0x6a8bf1e0);
}


TEST(rate_adaptive_code, no_ra_if_no_linecombs) {
    auto H = get_code_big_nora();
    EXPECT_ANY_THROW(H.set_rate(H.get_n_rows_mother_matrix() - 5));
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


// vn eliminations are allowed now! TODO reconsider this.
//TEST(rate_adaptive_code, dont_allow_vn_elimination) {
//    std::vector<std::uint32_t> colptr{0, 1, 2, 4, 5, 7, 9, 12};
//    std::vector<std::uint16_t> row_idx{0, 1, 0, 1, 2, 0, 2, 1, 2, 0, 1, 2};
//    EXPECT_ANY_THROW(RateAdaptiveCode<Bit>(colptr, row_idx, {0,1}));
//}


TEST(rate_adaptive_code, getters) {
    auto H = get_code_big_nora();
    EXPECT_EQ(H.get_n_rows_mother_matrix(), 2048);
    EXPECT_EQ(H.getNCols(), 6144);

//    H.get_current_n_rate_adapted_rows();
}

TEST(rate_adaptive_code, encode_with_ra) {
    auto H = get_code_big_wra();

    std::vector<Bit> input = get_bitstring(H.getNCols()); // true data to be sent

    // storage for syndrome. Initialize with arbitrary values, which must be overwritten by encoder.
    std::vector<Bit> syndrome = get_bitstring(H.get_n_rows_mother_matrix());

    // test agreement with encoder that doesn't use rate adaption
    H.encode_with_ra(input, syndrome, H.get_n_rows_mother_matrix());

    EXPECT_EQ(hash_vector(syndrome), 2814594723);

    // check that invalid requests lead to exceptions
    EXPECT_ANY_THROW(H.encode_with_ra({true, false}, syndrome, -1));  // invalid input size
    EXPECT_ANY_THROW(H.encode_with_ra(input, syndrome, -1)); // invalid requested size
    EXPECT_ANY_THROW(H.encode_with_ra(input, syndrome, H.get_n_rows_mother_matrix() + 1)); // too big requested size

    H.encode_with_ra(input, syndrome, H.get_n_rows_mother_matrix() / 2 + 1);
    EXPECT_EQ(hash_vector(syndrome), 0x4e395580);

    H.encode_with_ra(input, syndrome, static_cast<size_t>(H.get_n_rows_mother_matrix() * 0.7));
    EXPECT_EQ(hash_vector(syndrome), 0x01dab680);
}


TEST(rate_adaptive_code, ra_reported_size) {
    auto H = get_code_big_wra();

    {
        auto H_copy = H;
        H_copy.set_rate(0);
        EXPECT_EQ(H_copy, H);  // set_rate(0) does nothing.

        // rate adapting (5 steps) sets correct reported lengths
        constexpr std::size_t n_line_combs = 5;
        H_copy.set_rate(n_line_combs);
        EXPECT_EQ(H_copy.get_n_rows_after_rate_adaption(), H_copy.get_n_rows_mother_matrix() - n_line_combs);
    }
}

TEST(rate_adaptive_code, decode_infer_rate) {
    auto H = get_code_big_wra();

    std::vector<Bit> x = get_bitstring(H.getNCols()); // true data to be sent

    // storage for syndrome. Initialize with arbitrary values, which must be overwritten by encoder.
    std::vector<Bit> syndrome;
    constexpr double rate_adapt_factor = .95;
    H.encode_with_ra(x, syndrome, H.get_n_rows_mother_matrix() * rate_adapt_factor);

    constexpr double p = 0.005;
    std::vector<bool> x_noised = x; // copy for distorted data
    noise_bitstring_inplace(x_noised, p);
    ASSERT_FALSE(x_noised == x);  // actually have errors to be corrected!

    double vlog = log((1 - p) / p);
    std::vector<double> llrs(x.size());
    for (std::size_t i{}; i < llrs.size(); ++i) {
        llrs[i] = vlog * (1 - 2 * x_noised[i]); // log likelihood ratios
    }

    std::vector<bool> prediction;
    bool success = H.decode_infer_rate(llrs, syndrome, prediction);

    ASSERT_TRUE(success);
    ASSERT_EQ(prediction, x);
}


TEST(rate_adaptive_code, rate_adapted_fer) {
    // assert that the rate adapted FER (at set fraction of mother syndrome) is less than a certain value.
    // TODO performance seems pretty bad! Not consistent with AFF3CT results!
    std::mt19937_64 rng(42);
    auto H = get_code_big_wra();

    constexpr double p = 0.01;
    constexpr std::size_t num_frames_to_test = 100;
    constexpr std::uint16_t max_num_iter = 50;
    const auto syndrome_size = static_cast<std::size_t>(H.get_n_rows_mother_matrix() - 10);

    std::size_t num_frame_errors{};;
    std::size_t frame_idx{1};  // counts the number of iterations
    for (; frame_idx < num_frames_to_test; ++frame_idx) {
        std::vector<bool> x(H.getNCols()); // true data sent over a noisy channel
        noise_bitstring_inplace(rng, x, 0.5);  // choose it randomly.

        std::vector<bool> syndrome;  // syndrome for error correction, which is sent over a noise-less channel.
        H.encode_with_ra(x, syndrome, syndrome_size);

        std::vector<bool> x_noised = x; // copy for distorted data
        noise_bitstring_inplace(rng, x_noised, p);

        // log likelihood ratio (llr) computation
        double vlog = ::log((1 - p) / p);
        std::vector<double> llrs(x.size());
        for (std::size_t i{}; i < llrs.size(); ++i) {
            llrs[i] = vlog * (1 - 2 * x_noised[i]); // log likelihood ratios
        }

        std::vector<bool> solution;
        bool success = H.decode_infer_rate(llrs, syndrome, solution, max_num_iter);

        if (solution == x) {
            if (!success) {
                std::cerr << "DECODER GIVES CORRECT RESULT ALTHOUGH IT HAS NOT CONVERGED!!!!" << std::endl;
                FAIL();
            }
        } else {
            num_frame_errors++;
            if (success) {
                std::cerr << "\n\nDECODER CONVERGED TO WRONG CODEWORD!!!!\n" << std::endl;
                FAIL();
            }
        }
    }

    double fer = static_cast<double>(num_frame_errors) / static_cast<double>(num_frames_to_test);
    std::cout << "FER: " << fer << " ( " << num_frame_errors << " errors from " << num_frames_to_test << " frames )" << std::endl;
    ASSERT_TRUE(fer < 0.2);
}
