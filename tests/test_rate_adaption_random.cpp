//
// Created by Adomas Baliuka on 31.07.26.
//
//
// Unit tests for the random rate adaption generator (LDPC4QKD/rate_adaption_random.hpp).
//
// Golden values below were computed independently in Python (using `sympy` for primality/factorization,
// not this codebase), by directly porting the same algorithm description as `RateAdaptLCG::get_LCG_with_period`:
//   C = smallest value >= nextprime(m // 4) that is coprime to m;
//   A_minus_1 = radical(m); if m % 4 == 0 and A_minus_1 % 4 != 0: A_minus_1 *= 4; A = A_minus_1 + 1;
//   then iterate seed = (A * seed + C) % m.

#include <gtest/gtest.h>
#include <set>
#include <cstdint>
#include <concepts>
#include <limits>
#include <stdexcept>

#include "LDPC4QKD/rate_adaption_random.hpp"

using namespace LDPC4QKD;

namespace {
    std::vector<std::uint64_t> raw_lcg_sequence(std::uint64_t m, std::uint64_t seed, std::size_t count) {
        auto lcg = RateAdaptLCG::get_LCG_with_period(m, seed);
        std::vector<std::uint64_t> result(count);
        for (auto &v: result) {
            v = lcg.next();
        }
        return result;
    }

    /*!
     * Test-only utility (not part of the library's public API -- `rate_adaptive_code.hpp` consumes
     * `RateAdaptLCG::LCG::next()` directly, never a materialized array). Generates the full sequence of row-index pairs
     * an auto-generated (LCG-based) rate adaption would combine for a mother matrix with `n_mother_rows` rows, to
     * exercise `RateAdaptLCG::get_LCG_with_period` end-to-end and to build explicit-`rows_to_combine` codes for
     * comparison against auto-generated ones in `test_rate_adaptive_code.cpp`.
     *
     * The result is meant to be consumed with combined rows placed at the FRONT of the rate-adapted output (in the
     * order given here) and leftover rows placed after them in their original order -- this is the only placement
     * convention that allows the array to be consumed as a growing prefix as the number of combined pairs
     * increases (required for `RateAdaptiveCode::set_rate()`).
     *
     * @tparam idx_t unsigned integer type fitting the number of mother rows.
     * @param n_mother_rows number of rows in the mother parity check matrix.
     * @param seed seed for the underlying LCG. Defaults to 0, matching the reference implementation.
     * @return flat array of row-index pairs to combine: `result[2*i], result[2*i+1]` is the `i`-th pair.
     *         Has length `2 * (n_mother_rows / 2)`. Empty if `n_mother_rows < 2` (no combination possible).
     */
    template<std::unsigned_integral idx_t = std::uint16_t>
    std::vector<idx_t> generate_random_rate_adaption(std::size_t n_mother_rows, std::uint64_t seed = 0) {
        if (n_mother_rows < 2) {
            return {};
        }
        if (n_mother_rows - 1 > static_cast<std::size_t>(std::numeric_limits<idx_t>::max())) {
            throw std::domain_error(
                    "generate_random_rate_adaption: `idx_t` is too narrow to hold row indices for the given "
                    "`n_mother_rows` (mother row indices would silently truncate).");
        }

        const std::size_t max_line_combs = n_mother_rows / 2;
        auto lcg = RateAdaptLCG::get_LCG_with_period(n_mother_rows, seed);

        std::vector<idx_t> rows_to_combine(2 * max_line_combs);
        for (auto &v: rows_to_combine) {
            v = static_cast<idx_t>(lcg.next());
        }
        return rows_to_combine;
    }

    // Proves `get_LCG_with_period` is actually usable in a constant expression (not just marked `constexpr`
    // without it mattering): both the construction and `next()` are evaluated entirely at compile time here.
    constexpr std::uint64_t compute_first_value_at_compile_time(std::uint64_t m, std::uint64_t seed) {
        auto lcg = LDPC4QKD::RateAdaptLCG::get_LCG_with_period(m, seed);
        return lcg.next();
    }

    static_assert(compute_first_value_at_compile_time(2048, 0) == 521);
    static_assert(compute_first_value_at_compile_time(2, 0) == 1);
    static_assert(LDPC4QKD::RateAdaptLCG::get_LCG_with_period(6144, 0).A == 25);
}

// ---------------------------------------------------------------------------------------------- LCG golden values

TEST(rate_adaption_random, lcg_golden_values_m100_seed0) {
    std::vector<std::uint64_t> expect{29, 18, 67, 76, 45, 74, 63, 12};
    EXPECT_EQ(raw_lcg_sequence(100, 0, 8), expect);
}

TEST(rate_adaption_random, lcg_golden_values_m2048_seed0) {
    std::vector<std::uint64_t> expect{521, 1114, 307, 1236, 1405, 878, 231, 552};
    EXPECT_EQ(raw_lcg_sequence(2048, 0, 8), expect);
}

TEST(rate_adaption_random, lcg_golden_values_m2048_seed42) {
    // Confirms the `seed` parameter is actually threaded through (not just always starting at 0).
    std::vector<std::uint64_t> expect{899, 420, 205, 318, 1335, 248, 705, 722};
    EXPECT_EQ(raw_lcg_sequence(2048, 42, 8), expect);
}

TEST(rate_adaption_random, lcg_golden_values_m6144_seed0) {
    std::vector<std::uint64_t> expect{1543, 3254, 3021, 3340, 5171, 1794, 3385, 152};
    EXPECT_EQ(raw_lcg_sequence(6144, 0, 8), expect);
}

TEST(rate_adaption_random, lcg_golden_values_odd_m6143_seed0) {
    // Sanity check that odd `m` (not used in practice for mother matrices, but not disallowed either) works.
    std::vector<std::uint64_t> expect{1543, 3086, 4629, 29, 1572, 3115, 4658, 58};
    EXPECT_EQ(raw_lcg_sequence(6143, 0, 8), expect);
}

// ------------------------------------------------------ full period for every `m` (no more degenerate `m` cases)

// m=2, m=4, m=10 used to be exactly the cases where naively using `C = nextprime(m/4)` fails to be coprime to
// `m` (e.g. nextprime(0) = 2, and gcd(2, 2) = 2). `get_LCG_with_period` now searches forward from that
// candidate for one that actually is coprime to `m`, so these succeed and give a genuinely full-period
// sequence instead of throwing. Golden values cross-checked independently in Python.

TEST(rate_adaption_random, lcg_full_period_m2) {
    std::vector<std::uint64_t> expect{1, 0};
    const auto seq = raw_lcg_sequence(2, 0, 2);
    EXPECT_EQ(seq, expect);
    EXPECT_EQ(std::set<std::uint64_t>(seq.begin(), seq.end()).size(), 2u);
}

TEST(rate_adaption_random, lcg_full_period_m4) {
    std::vector<std::uint64_t> expect{3, 2, 1, 0};
    const auto seq = raw_lcg_sequence(4, 0, 4);
    EXPECT_EQ(seq, expect);
    EXPECT_EQ(std::set<std::uint64_t>(seq.begin(), seq.end()).size(), 4u);
}

TEST(rate_adaption_random, lcg_full_period_m10) {
    std::vector<std::uint64_t> expect{3, 6, 9, 2, 5, 8, 1, 4, 7, 0};
    const auto seq = raw_lcg_sequence(10, 0, 10);
    EXPECT_EQ(seq, expect);
    EXPECT_EQ(std::set<std::uint64_t>(seq.begin(), seq.end()).size(), 10u);
}

TEST(rate_adaption_random, lcg_full_period_holds_generally) {
    // Beyond the three specific golden cases above: for a range of `m` (including ones that were never
    // degenerate to begin with), confirm a full cycle of `m` draws visits every residue in [0, m) exactly once.
    for (const std::uint64_t m: {2ul, 3ul, 4ul, 5ul, 6ul, 7ul, 8ul, 9ul, 10ul, 11ul, 12ul, 100ul, 2048ul, 6144ul}) {
        const auto seq = raw_lcg_sequence(m, 0, m);
        const std::set<std::uint64_t> unique_vals(seq.begin(), seq.end());
        EXPECT_EQ(unique_vals.size(), m) << "m=" << m;
        EXPECT_EQ(*unique_vals.begin(), 0u) << "m=" << m;
        EXPECT_EQ(*unique_vals.rbegin(), m - 1) << "m=" << m;
    }
}

TEST(rate_adaption_random, lcg_m1_is_trivial) {
    // Degenerate but well-defined: the only residue modulo 1 is 0.
    auto lcg = RateAdaptLCG::get_LCG_with_period(1, 0);
    EXPECT_EQ(lcg.next(), 0u);
    EXPECT_EQ(lcg.next(), 0u);
}

TEST(rate_adaption_random, lcg_zero_m_throws) {
    EXPECT_THROW(RateAdaptLCG::get_LCG_with_period(0, 0), std::domain_error);
}

// --------------------------------------------------------------------------------- generate_random_rate_adaption

TEST(rate_adaption_random, small_n_mother_rows_returns_empty) {
    EXPECT_TRUE(generate_random_rate_adaption<std::uint16_t>(0).empty());
    EXPECT_TRUE(generate_random_rate_adaption<std::uint16_t>(1).empty());
}

TEST(rate_adaption_random, idx_t_too_narrow_throws) {
    // n_mother_rows=1000 does not fit in uint8_t (max 255): must throw, not silently truncate row indices.
    EXPECT_THROW(generate_random_rate_adaption<std::uint8_t>(1000), std::domain_error);
    // but a value that does fit should not throw:
    EXPECT_NO_THROW(generate_random_rate_adaption<std::uint8_t>(200));
}

TEST(rate_adaption_random, length_matches_max_line_combs) {
    for (const std::size_t m: {2048ul, 4096ul, 6144ul, 8192ul, 16384ul, 24576ul}) {
        const auto rows_to_combine = generate_random_rate_adaption<std::uint16_t>(m);
        EXPECT_EQ(rows_to_combine.size(), 2 * (m / 2)) << "m=" << m;
    }
    // odd m: one row is never used (floor division truncates the pairing).
    const auto rows_to_combine_odd = generate_random_rate_adaption<std::uint16_t>(6143);
    EXPECT_EQ(rows_to_combine_odd.size(), 2 * (6143 / 2));
}

TEST(rate_adaption_random, all_indices_within_bounds) {
    for (const std::size_t m: {2048ul, 6144ul, 6143ul, 24576ul}) {
        const auto rows_to_combine = generate_random_rate_adaption<std::uint16_t>(m);
        for (const auto idx: rows_to_combine) {
            EXPECT_LT(idx, m) << "m=" << m;
        }
    }
}

TEST(rate_adaption_random, no_duplicate_indices) {
    // A full-period LCG visits every value in [0, m) at most once over m steps, so a prefix of `2*(m/2) <= m`
    // draws must consist of distinct values (no mother row is combined twice).
    for (const std::size_t m: {2048ul, 6144ul, 6143ul, 24576ul}) {
        const auto rows_to_combine = generate_random_rate_adaption<std::uint16_t>(m);
        const std::set<std::uint16_t> unique_indices(rows_to_combine.begin(), rows_to_combine.end());
        EXPECT_EQ(unique_indices.size(), rows_to_combine.size()) << "m=" << m;
    }
}

TEST(rate_adaption_random, matches_raw_lcg_sequence) {
    // `generate_random_rate_adaption` should be nothing more than "pull 2*(m/2) values from the same LCG",
    // in order -- this pins down that relationship so the two can't silently drift apart.
    constexpr std::size_t m = 6144;
    const auto rows_to_combine = generate_random_rate_adaption<std::uint16_t>(m);
    const auto raw = raw_lcg_sequence(m, 0, 2 * (m / 2));
    ASSERT_EQ(rows_to_combine.size(), raw.size());
    for (std::size_t i = 0; i < raw.size(); ++i) {
        EXPECT_EQ(rows_to_combine[i], raw[i]) << "i=" << i;
    }
}

TEST(rate_adaption_random, determinism_same_seed_same_result) {
    const auto a = generate_random_rate_adaption<std::uint16_t>(2048, 7);
    const auto b = generate_random_rate_adaption<std::uint16_t>(2048, 7);
    EXPECT_EQ(a, b);
}

TEST(rate_adaption_random, different_seeds_give_different_results) {
    const auto a = generate_random_rate_adaption<std::uint16_t>(2048, 0);
    const auto b = generate_random_rate_adaption<std::uint16_t>(2048, 42);
    EXPECT_NE(a, b);
}

// ------------------------------------------------------------------- real mother-matrix sizes used in this repo

TEST(rate_adaption_random, real_mother_sizes_do_not_hit_degenerate_case) {
    // These are the actual mother ROW counts (= M * expansion_factor, not the column count N despite the
    // "MxN"-style naming of the encoder namespaces) of every code in this repository:
    // - protograph-based: 2048 (`get_rate_adaptive_code(0)` and `(3)`), 8192 (`(1)` and `(4)`),
    //   524288 (`(2)` and `(5)`).
    // - 819k degree-distribution-based: lrate 0.1 through 0.5, i.e. `get_code_819k(6..14)`.
    // None of them may hit the degenerate (non-coprime `C`) edge case, or the fallback rate adaption used by
    // `RateAdaptiveCode` would throw at construction time for these real codes.
    for (const std::size_t m: {2048ul, 8192ul, 524288ul,
                                81920ul, 122880ul, 163840ul, 204800ul, 245760ul, 286720ul, 327680ul, 368640ul,
                                409600ul}) {
        EXPECT_NO_THROW(generate_random_rate_adaption<std::uint32_t>(m)) << "m=" << m;
    }
}

TEST(rate_adaption_random, golden_values_m204800_seed0) {
    // m = 204800 is the mother row count of `lrate_0.25_block_819k` (`get_code_819k(9)`). These values, combined
    // with the front-placed union-of-mother-rows scheme in `rate_adaptive_code.hpp`, were confirmed to reproduce
    // the externally/independently-generated `alist_03ea702_ra_Hra1.alist` and `..._Hra2.alist` fixtures
    // bit-for-bit, up to row order (both files live in `TMPTMP/`, see `rate_adaption_plan.md` for the scheme
    // description; the fixtures identify their source mother matrix by content hash, not by this repo's naming).
    std::vector<std::uint64_t> expect{51203, 102526, 158769, 7132, 138815, 8218, 183341, 195384};
    EXPECT_EQ(raw_lcg_sequence(204800, 0, 8), expect);
}

TEST(rate_adaption_random, golden_values_m81920_seed0) {
    // m = 81920 is the mother row count of `lrate_0.1_block_819k` (`get_code_819k(6)`). Same validation as
    // `golden_values_m204800_seed0` above, against `alist_d1f3a91_ra_Hra1.alist` / `..._Hra2.alist`.
    std::vector<std::uint64_t> expect{20483, 41086, 66609, 48092, 26175, 28698, 50221, 31544};
    EXPECT_EQ(raw_lcg_sequence(81920, 0, 8), expect);
}
