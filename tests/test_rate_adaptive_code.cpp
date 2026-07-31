//
// Created by alice on, 07.05.21.
//

// Google Test framework
#include <gtest/gtest.h>
#include "helpers_for_testing.hpp"

// Standard library
#include <iostream>

// To be tested
#include "LDPC4QKD/rate_adaptive_code.hpp"
#include "LDPC4QKD/prebuilt_codes.hpp"
#include "LDPC4QKD/rate_adaption_random.hpp"

// Test cases test against constants known to be correct for the LDPC-matrix defined here:
#include "fortest_autogen_ldpc_matrix_csc.hpp"
#include "fortest_autogen_rate_adaption.hpp"

using namespace HelpersForTests;
using namespace LDPC4QKD;

namespace {

    auto get_code_big_nora() {
        std::vector<std::uint32_t> colptr(AutogenLDPC::colptr.begin(), AutogenLDPC::colptr.end());
        std::vector<std::uint16_t> row_idx(AutogenLDPC::row_idx.begin(), AutogenLDPC::row_idx.end());
        return RateAdaptiveCode(colptr, row_idx);
    }

    auto get_code_big_wra() {
        std::vector<std::uint32_t> colptr(AutogenLDPC::colptr.begin(), AutogenLDPC::colptr.end());
        std::vector<std::uint16_t> row_idx(AutogenLDPC::row_idx.begin(), AutogenLDPC::row_idx.end());
        std::vector<std::uint16_t> rows_to_combine(AutogenRateAdapt::rows.begin(), AutogenRateAdapt::rows.end());
        return RateAdaptiveCode(colptr, row_idx, rows_to_combine);
    }

    auto get_code_small() {
        //    H =  [1 0 1 0 1 0 1
        //			0 1 1 0 0 1 1
        //			0 0 0 1 1 1 1]
        std::vector<std::uint32_t> colptr{0, 1, 2, 4, 5, 7, 9, 12};
        std::vector<std::uint16_t> row_idx{0, 1, 0, 1, 2, 0, 2, 1, 2, 0, 1, 2};
        return RateAdaptiveCode(colptr, row_idx);
    }

    auto get_code_819k(std::size_t id) {
        return HelperFixedSize::get_rate_adaptive_code(id);
    }

}

TEST(rate_adaptive_code_from_colptr_rowIdx, decode_test_small) {
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

TEST(rate_adaptive_code_from_colptr_rowIdx, decode_test_big) {
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


TEST(rate_adaptive_code_from_colptr_rowIdx, encode_no_ra) {
    auto H = get_code_big_nora();
    std::vector<Bit> in = get_bitstring(H.getNCols());
    std::vector<Bit> out(H.get_n_rows_mother_matrix());

    std::cout << hash_vector(in) << std::endl;
    H.encode_no_ra(in, out);

    EXPECT_EQ(hash_vector(out), 2814594723);
}


TEST(rate_adaptive_code_from_colptr_rowIdx, encode_current_rate) {
    auto H = get_code_big_wra();
    std::vector<Bit> in = get_bitstring(H.getNCols());
    std::vector<Bit> out(H.get_n_rows_mother_matrix());

    H.encode_at_current_rate(in, out);

    EXPECT_EQ(hash_vector(out), 2814594723);

    const auto n_line_combs = static_cast<std::size_t>(static_cast<double>(H.get_n_rows_mother_matrix()) * 0.3);
    H.set_rate(n_line_combs);

    H.encode_at_current_rate(in, out);
    EXPECT_EQ(hash_vector(out), 0x6a8bf1e0);
}


TEST(rate_adaptive_code_from_colptr_rowIdx, set_rate_throws_beyond_max_steps) {
    // mother-only-constructed code: auto-generated rate adaption supports up to `n_mother_rows / 2` line combs.
    auto H = get_code_big_nora();
    EXPECT_ANY_THROW(H.set_rate(H.get_max_ra_steps() + 1));
}

TEST(rate_adaptive_code_from_colptr_rowIdx, set_rate_at_max_steps_succeeds) {
    auto H = get_code_big_nora();
    EXPECT_NO_THROW(H.set_rate(H.get_max_ra_steps()));
    EXPECT_EQ(H.get_n_rows_after_rate_adaption(), H.get_n_rows_mother_matrix() - H.get_max_ra_steps());
}

TEST(rate_adaptive_code_from_colptr_rowIdx, auto_rate_adaption_max_steps) {
    // Before this change, a mother-only-constructed code had `rows_to_combine` empty AND no rate adaption at all
    // (`get_max_ra_steps() == 0`). Now it has real auto-generated (LCG-based) rate adaption instead.
    auto H = get_code_big_nora();
    EXPECT_EQ(H.get_max_ra_steps(), H.get_n_rows_mother_matrix() / 2);
}


TEST(rate_adaptive_code_from_colptr_rowIdx, init_pos_CN_pos_VN) {
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


TEST(rate_adaptive_code_from_colptr_rowIdx, auto_rate_adaption_golden_small) {
    // `get_code_small()` is mother-only-constructed (no explicit `rows_to_combine`), so `set_rate` must use the
    // auto-generated (LCG-based) scheme. n_mother_rows = 3, so LCG(m=3, seed=0) applies; golden values for this
    // LCG were independently cross-checked in Python (see test_rate_adaption_random.cpp for the same method):
    // it yields the pair (2, 1) -- i.e. combine mother rows 2 and 1.
    auto H = get_code_small();
    ASSERT_EQ(H.get_max_ra_steps(), 1u); // n_mother_rows / 2 = 3 / 2 = 1

    H.set_rate(1);

    // front-placed convention: the combined row goes to position 0, the leftover row (0) to position 1.
    // mother row 1 = {1, 2, 5, 6}, mother row 2 = {3, 4, 5, 6} (see `init_pos_CN_pos_VN` above);
    // their union (sorted, deduplicated) is {1, 2, 3, 4, 5, 6}. Mother row 0 = {0, 2, 4, 6} is the leftover.
    std::vector<std::vector<decltype(H)::MatrixIndex>> expect_posVN{{1, 2, 3, 4, 5, 6},
                                                                    {0, 2, 4, 6}};
    EXPECT_EQ(H.getPosVarn(), expect_posVN);
}


TEST(rate_adaptive_code_from_colptr_rowIdx, auto_rate_adaption_round_trip) {
    // mother-only-constructed code (no explicit `rows_to_combine`): verify encode/decode agree across a range of
    // auto-generated (LCG-based) rate-adapted rates.
    auto H = get_code_big_nora();
    std::vector<Bit> x = get_bitstring(H.getNCols());

    constexpr double p = 0.02;
    std::vector<Bit> x_noised = x;
    noise_bitstring_inplace(x_noised, p);
    ASSERT_FALSE(x_noised == x);

    double vlog = log((1 - p) / p);
    std::vector<double> llrs(x.size());
    for (std::size_t i{}; i < llrs.size(); ++i) {
        llrs[i] = vlog * (1 - 2 * x_noised[i]);
    }

    for (const double frac: {1.0, 0.95, 0.9, 0.8, 0.7}) {
        const auto syndrome_len = static_cast<std::size_t>(
                static_cast<double>(H.get_n_rows_mother_matrix()) * frac);

        std::vector<Bit> syndrome;
        H.encode_with_ra(x, syndrome, syndrome_len);

        std::vector<Bit> solution;
        const bool success = H.decode_infer_rate(llrs, syndrome, solution);
        EXPECT_TRUE(success) << "frac=" << frac;
        EXPECT_EQ(solution, x) << "frac=" << frac;
    }
}


TEST(rate_adaptive_code_from_colptr_rowIdx, auto_rate_adaption_equals) {
    auto H1 = get_code_big_nora();
    auto H2 = get_code_big_nora();
    EXPECT_TRUE(H1 == H2); // two auto-generated codes from the same mother compare equal.
    H1.set_rate(3);
    H2.set_rate(3);
    EXPECT_TRUE(H1 == H2);

    // explicitly feeding back the auto-generated scheme's own materialized pairs is NOT equivalent: a non-empty
    // `rows_to_combine` always uses the back-placed convention (see rate_adaption_random.hpp caveat). Materialize
    // those pairs by hand here (rather than via a library-level utility -- there is none; `rate_adaptive_code.hpp`
    // consumes `RateAdaptLCG::LCG::next()` directly, never a materialized array) by replaying the same LCG the
    // auto-generated scheme uses internally.
    std::vector<std::uint32_t> colptr(AutogenLDPC::colptr.begin(), AutogenLDPC::colptr.end());
    std::vector<std::uint16_t> row_idx(AutogenLDPC::row_idx.begin(), AutogenLDPC::row_idx.end());
    auto lcg = RateAdaptLCG::get_LCG_with_period(H1.get_n_rows_mother_matrix(), 0);
    std::vector<std::uint16_t> explicit_pairs(2 * (H1.get_n_rows_mother_matrix() / 2));
    for (auto &v: explicit_pairs) {
        v = static_cast<std::uint16_t>(lcg.next());
    }
    RateAdaptiveCode<std::uint16_t> H3(colptr, row_idx, explicit_pairs);

    EXPECT_FALSE(H1 == H3);
}


TEST(rate_adaptive_code_from_decoder, auto_rate_adaption_no_materialized_array) {
    // A mother-only-constructed 819k-sized code must not materialize an O(n_mother_rows)-sized `rows_to_combine`
    // vector. There is no direct public getter for `rows_to_combine`, but `get_max_ra_steps()` can only equal
    // `n_mother_rows / 2` via the auto-generated (empty `rows_to_combine`) branch: `prebuilt_codes.hpp` passes an
    // explicit empty `std::vector<Idx>{}` for these codes (see `get_rate_adaptive_code` cases 6-14), so this
    // indirectly confirms no large vector was built.
    auto H = get_code_819k(6);
    EXPECT_EQ(H.get_max_ra_steps(), H.get_n_rows_mother_matrix() / 2);
}


TEST(rate_adaptive_code_from_decoder, get_code_819k_round_trip_with_rate_adaption) {
    // Cases 6-13 (819k degree-distribution codes) previously passed an explicit empty `rows_to_combine`, which
    // used to mean NO rate adaption at all. Now this triggers real auto-generated (LCG-based) rate adaption;
    // verify encode/decode still agree once some lines are actually combined (unlike `new_819k_code.fer_simulation`,
    // which only tests case 14 at the mother rate, without exercising rate adaption).
    // Only the smallest (6) and largest (13) mother matrices are exercised here: the rate-adaption code path under
    // test (front-placed LCG scheme in `recompute_pos_vn_cn`/`encode_with_ra`) is identical for every id, and
    // doesn't depend on matrix content -- only on `n_mother_rows`, whose generality across all 9 819k mother row
    // counts is separately (and much more cheaply) covered by
    // `rate_adaption_random.real_mother_sizes_do_not_hit_degenerate_case`. Testing all 8 ids here would mostly
    // re-run the same code path at extra BP-decode cost without adding independent coverage.
    // `p` must stay well under the channel capacity bound of the least redundant code tested here (case 6,
    // lrate 0.1: capacity bound is at a bit-flip probability of ~0.013), which is far more restrictive than the
    // 0.03-0.5 range used for the other (much higher-lrate) codes tested elsewhere in this file.
    std::mt19937_64 rng(42);
    constexpr double p = 0.006;
    constexpr std::size_t max_num_iter = 50;

    for (const std::size_t id: {6u, 13u}) {
        auto H = get_code_819k(id);
        const std::size_t syndrome_size = H.get_n_rows_mother_matrix() - 10;

        std::vector<bool> x(H.getNCols());
        noise_bitstring_inplace(rng, x, 0.5);

        std::vector<bool> syndrome;
        H.encode_with_ra(x, syndrome, syndrome_size);

        std::vector<bool> x_noised = x;
        noise_bitstring_inplace(rng, x_noised, p);

        double vlog = ::log((1 - p) / p);
        std::vector<double> llrs(x.size());
        for (std::size_t i{}; i < llrs.size(); ++i) {
            llrs[i] = vlog * (1 - 2 * x_noised[i]);
        }

        std::vector<bool> solution;
        const bool success = H.decode_infer_rate(llrs, syndrome, solution, max_num_iter);
        EXPECT_TRUE(success) << "id=" << id;
        EXPECT_EQ(solution, x) << "id=" << id;
    }
}


// vn eliminations are allowed now! TODO reconsider this.
//TEST(rate_adaptive_code_from_colptr_rowIdx, dont_allow_vn_elimination) {
//    std::vector<std::uint32_t> colptr{0, 1, 2, 4, 5, 7, 9, 12};
//    std::vector<std::uint16_t> row_idx{0, 1, 0, 1, 2, 0, 2, 1, 2, 0, 1, 2};
//    EXPECT_ANY_THROW(RateAdaptiveCode(colptr, row_idx, {0,1}));
//}


TEST(rate_adaptive_code_from_colptr_rowIdx, getters) {
    auto H = get_code_big_nora();
    EXPECT_EQ(H.get_n_rows_mother_matrix(), 2048);
    EXPECT_EQ(H.getNCols(), 6144);

    EXPECT_EQ(H.get_n_rows_after_rate_adaption(), H.get_n_rows_mother_matrix());
}

TEST(rate_adaptive_code_from_colptr_rowIdx, encode_with_ra) {
    auto H = get_code_big_wra();

    std::vector<Bit> input = get_bitstring(H.getNCols()); // true data to be sent

    // storage for syndrome. Initialize with arbitrary values, which must be overwritten by encoder.
    std::vector<Bit> syndrome = get_bitstring(H.get_n_rows_mother_matrix());

    // test agreement with encoder that doesn't use rate adaption
    H.encode_with_ra(input, syndrome, H.get_n_rows_mother_matrix());

    EXPECT_EQ(hash_vector(syndrome), 2814594723);

    // check that invalid requests lead to exceptions
    EXPECT_ANY_THROW(H.encode_with_ra({true, false}, syndrome, 0xbadb'eeff'ffff'aaaa));  // invalid input size
    EXPECT_ANY_THROW(H.encode_with_ra(input, syndrome, 0xbadb'eeff'ffff'ffff)); // invalid requested size
    EXPECT_ANY_THROW(H.encode_with_ra(input, syndrome, H.get_n_rows_mother_matrix() + 1)); // too big requested size

    H.encode_with_ra(input, syndrome, H.get_n_rows_mother_matrix() / 2 + 1);
    EXPECT_EQ(hash_vector(syndrome), 0x4e395580);

    H.encode_with_ra(input, syndrome,
                     static_cast<size_t>(static_cast<double>(H.get_n_rows_mother_matrix()) * 0.7));
    EXPECT_EQ(hash_vector(syndrome), 0x01dab680);
}


TEST(rate_adaptive_code_from_colptr_rowIdx, ra_reported_size) {
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

TEST(rate_adaptive_code_from_colptr_rowIdx, decode_infer_rate) {
    auto H = get_code_big_wra();

    std::vector<Bit> x = get_bitstring(H.getNCols()); // true data to be sent

    // storage for syndrome. Initialize with arbitrary values, which must be overwritten by encoder.
    std::vector<Bit> syndrome;
    constexpr double rate_adapt_factor = .95;
    const auto output_syndr_len = static_cast<size_t>(
            static_cast<double>(H.get_n_rows_mother_matrix()) * rate_adapt_factor);
    H.encode_with_ra(x, syndrome, output_syndr_len);

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

TEST(rate_adaptive_code_819k, fer_simulation) {
    // assert that the rate adapted FER (at set fraction of mother syndrome) is small.
    std::mt19937_64 rng(42);
    auto H = get_code_819k(14);

    constexpr double p = 0.095;
    constexpr std::size_t num_frames_to_test = 2;
    constexpr std::size_t max_num_iter = 50;
    const std::size_t syndrome_size = H.get_n_rows_mother_matrix();

    std::cout << "Testing code with size " << H.get_n_rows_mother_matrix() << " x " << H.getNCols()
              << " at QBER = " << p << std::endl;

    std::size_t num_frame_errors{};
    std::size_t frame_idx{0};  // counts the number of iterations
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
    std::cout << "FER: " << fer << " ( " << num_frame_errors << " errors from " << num_frames_to_test << " frames )"
              << std::endl;
    ASSERT_EQ(fer, 0.);
}

TEST(rate_adaptive_code_from_colptr_rowIdx, rate_adapted_fer) {
    // assert that the rate adapted FER (at set fraction of mother syndrome) is small.
    std::mt19937_64 rng(42);
    auto H = get_code_big_wra();

    constexpr double p = 0.03;
    constexpr std::size_t num_frames_to_test = 100;
    constexpr std::size_t max_num_iter = 50;
    const std::size_t syndrome_size = H.get_n_rows_mother_matrix() - 10;

    std::cout << "Testing rate adaption for mother code with size " << H.get_n_rows_mother_matrix() << " x " << H.getNCols()
              << " at QBER = " << p << std::endl;

    std::size_t num_frame_errors{};
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
    std::cout << "FER: " << fer << " ( " << num_frame_errors << " errors from " << num_frames_to_test << " frames )"
              << std::endl;
    ASSERT_EQ(fer, 0.);
}

TEST(rate_adaptive_code_from_colptr_rowIdx, llrs_bsc) {
    std::vector<bool> x{1, 1, 1, 1, 0, 0, 0};
    double p = 0.01;

    double vlog = log((1 - p) / p);
    std::vector<double> llrs(x.size());
    for (std::size_t i{}; i < llrs.size(); ++i) {
        llrs[i] = vlog * (1 - 2 * x[i]); // log likelihood ratios
    }

    std::vector<double> llrs_convenience = LDPC4QKD::llrs_bsc(x, p);

    EXPECT_EQ(llrs_convenience, llrs);
}

TEST(rate_adaptive_code_from_colptr_rowIdx, equals_not_equals_operators) {
    auto H1 = get_code_big_wra();
    auto H2 = get_code_big_wra();
    EXPECT_FALSE(H1 != H2);
    EXPECT_TRUE(H1 == H2);
    H1.set_rate(1);
    EXPECT_FALSE(H1 == H2);
    EXPECT_TRUE(H1 != H2);
}

TEST(rate_adaptive_code_from_decoder, obtain_from_advanced_encoder_behaviour) {
    std::vector<std::uint16_t> rows_to_combine{}; // not used here!
    RateAdaptiveCode<std::uint16_t> H1(encoder_2048x6144_4663d91.get_pos_varn(), rows_to_combine);

    auto H2 = get_code_big_wra();

    {
        std::vector<std::uint8_t> in = get_bitstring<std::uint8_t>(H1.getNCols());
        std::vector<std::uint8_t> out(H1.get_n_rows_mother_matrix());

        std::cout << "input hash: " << hash_vector(in) << std::endl;
        H1.encode_no_ra(in, out);
        std::cout << "output hash: " << hash_vector(out) << std::endl;

        EXPECT_EQ(hash_vector(out), 2814594723);
    }
    {
        std::vector<std::uint8_t> in = get_bitstring<std::uint8_t>(H2.getNCols());
        std::vector<std::uint8_t> out(H2.get_n_rows_mother_matrix());

        std::cout << "input hash: " << hash_vector(in) << std::endl;
        H2.encode_no_ra(in, out);
        std::cout << "output hash: " << hash_vector(out) << std::endl;

        EXPECT_EQ(hash_vector(out), 2814594723);
    }
}

TEST(rate_adaptive_code_from_decoder, obtain_from_advanced_encoder_equals) {
    std::vector<std::uint16_t> rows_to_combine(AutogenRateAdapt::rows.begin(), AutogenRateAdapt::rows.end());
    RateAdaptiveCode<std::uint16_t> H1(encoder_2048x6144_4663d91.get_pos_varn(), rows_to_combine);

    // TODO add random rate adaption for comparison
    auto H2 = get_code_big_wra();
    EXPECT_TRUE(H1 == H2);
}
