//
// Tests for the generic, compile-time QC-LDPC encoder machinery in `fixed_size_encoder.hpp`.
// For tests of the concrete prebuilt codes built on top of this machinery, see `test_prebuilt_codes.cpp`.
//

#include <gtest/gtest.h>

#include <array>
#include <vector>

#include "LDPC4QKD/fixed_size_encoder.hpp"
#include "LDPC4QKD/prebuilt_codes.hpp"
#include "LDPC4QKD/rate_adaptive_code.hpp"

// Ground truth for the same 2048x6144 matrix as `encoder_2048x6144_4663d91`, stored in full
// (expanded) CSC form rather than QC-exponent form.
#include "fortest_autogen_ldpc_matrix_csc.hpp"

using namespace LDPC4QKD;

// `bits_needed`, `smallest_type`, and `xor_as_bools` are all `constexpr`, so these are checked at
// compile time rather than as gtest runtime tests. `FixedSizeEncoder`/`FixedSizeEncoderQC`'s member
// functions are all `constexpr` too, so the toy-matrix checks further below are also static_asserts.
// The two remaining runtime TESTs check thrown exceptions (inherently runtime-only: a throw during
// constant evaluation is a compile error, not something a static_assert can observe) and cross-check
// against `RateAdaptiveCode`, which isn't `constexpr`-ready.
static_assert(bits_needed<0>() == 0);
static_assert(bits_needed<1>() == 1);
static_assert(bits_needed<2>() == 2);
static_assert(bits_needed<255>() == 8);
static_assert(bits_needed<256>() == 9);
static_assert(bits_needed<65535>() == 16);
static_assert(bits_needed<65536>() == 17);

static_assert(std::is_same_v<smallest_type<255>, std::uint8_t>);
static_assert(std::is_same_v<smallest_type<256>, std::uint16_t>);
static_assert(std::is_same_v<smallest_type<65535>, std::uint16_t>);
static_assert(std::is_same_v<smallest_type<65536>, std::uint32_t>);

static_assert(!xor_as_bools(0, 0));
static_assert(xor_as_bools(0, 1));
static_assert(xor_as_bools(1, 0));
static_assert(!xor_as_bools(1, 1));
// non-zero values other than 1 are also treated as "true"
static_assert(xor_as_bools(5, 0));
static_assert(!xor_as_bools(5, 3));

namespace {
    // A single 4x4 circulant block (a 1x1 matrix of exponents, expansion_factor=4), shifted by 2.
    // This is a permutation matrix equivalent to a cyclic shift of the input by 2 positions:
    // in[0]->out[2], in[1]->out[3], in[2]->out[0], in[3]->out[1].
    using ToyEncoder = FixedSizeEncoderQC<std::uint8_t, /*M=*/1, /*N=*/1, /*expansion_factor=*/4, /*num_nz=*/1,
            std::uint8_t, std::uint8_t, std::uint8_t>;

    constexpr ToyEncoder make_toy_encoder() {
        return ToyEncoder{
                std::array<std::uint8_t, 2>{0, 1},  // colptr
                std::array<std::uint8_t, 1>{0},     // row_idx
                std::array<std::uint8_t, 1>{2},     // values (shift)
        };
    }

    constexpr bool check_encode_qc() {
        auto enc = make_toy_encoder();
        std::array<std::uint8_t, 4> in{1, 0, 1, 1};
        std::array<std::uint8_t, 4> out{};
        enc.encode_qc(std::span<std::uint8_t const, 4>{in}, std::span<std::uint8_t, 4>{out});
        return out == std::array<std::uint8_t, 4>{1, 1, 1, 0};
    }

    constexpr bool check_get_pos_varn() {
        auto enc = make_toy_encoder();
        auto pos_varn = enc.get_pos_varn();
        return pos_varn.size() == 4
               && pos_varn[2] == std::vector<decltype(enc)::idx_t>{0}
               && pos_varn[3] == std::vector<decltype(enc)::idx_t>{1}
               && pos_varn[0] == std::vector<decltype(enc)::idx_t>{2}
               && pos_varn[1] == std::vector<decltype(enc)::idx_t>{3};
    }
}

static_assert(make_toy_encoder().get_input_size() == 4);
static_assert(make_toy_encoder().get_output_size() == 4);
static_assert(check_encode_qc());
static_assert(check_get_pos_varn());

// GCC currently can't constant-evaluate this specific overload (the generic, runtime-size-checked
// `encode()` called with `std::vector` arguments) even though it's declared `constexpr` -- it hits
// a compiler limitation around the virtual `encode_span` call from within this templated base-class
// member function. `check_encode_qc()` above already exercises the same underlying encoding logic
// via the non-throwing `encode_qc(span, span)` overload at compile time, so this stays a runtime
// test mainly to cover the `std::vector` / runtime-size-check code path itself.
TEST(test_fixed_size_encoder, qc_encoder_encode_generic_runtime_checked) {
    auto enc = make_toy_encoder();
    std::vector<std::uint8_t> in{1, 0, 1, 1};
    std::vector<std::uint8_t> out(4);
    enc.encode(in, out);
    EXPECT_EQ(out, (std::vector<std::uint8_t>{1, 1, 1, 0}));
}

TEST(test_fixed_size_encoder, qc_encoder_wrong_size_throws) {
    auto enc = make_toy_encoder();
    std::vector<std::uint8_t> in{1, 0, 1};  // wrong size: should be 4
    std::vector<std::uint8_t> out(4);
    EXPECT_THROW(enc.encode(in, out), std::out_of_range);
}

TEST(test_fixed_size_encoder, qc_encoder_matches_ground_truth_csc_matrix) {
    // Cross-validates `FixedSizeEncoderQC::get_pos_varn()` against a real, non-trivial matrix
    // (multiple nonzeros per column, real QC block structure) rather than the toy example above.
    // `encoder_2048x6144_4663d91` (QC-exponent form) and `AutogenLDPC` (full CSC form) are known to
    // represent the exact same LDPC matrix (see `test_rate_adaptive_code.cpp`'s
    // `obtain_from_advanced_encoder_behaviour`, which shows both give the same encode-hash).
    std::vector<std::uint32_t> colptr(AutogenLDPC::colptr.begin(), AutogenLDPC::colptr.end());
    std::vector<std::uint16_t> row_idx(AutogenLDPC::row_idx.begin(), AutogenLDPC::row_idx.end());
    RateAdaptiveCode<std::uint16_t> ground_truth(colptr, row_idx);

    const auto &expected_pos_varn = ground_truth.getPosVarn();
    auto actual_pos_varn = encoder_2048x6144_4663d91.get_pos_varn();

    ASSERT_EQ(actual_pos_varn.size(), expected_pos_varn.size());
    for (std::size_t row = 0; row < expected_pos_varn.size(); ++row) {
        EXPECT_EQ(actual_pos_varn[row], expected_pos_varn[row]) << "mismatch at row " << row;
    }
}

TEST(test_fixed_size_encoder, qc_encoder_invalid_matrix_throws_at_runtime) {
    // row_idx references row 99 of the matrix of exponents, but M=1 means that matrix only has
    // row 0: this must fail the out-of-memory-access check in the constructor.
    // Constructed as a plain (non-constexpr) runtime object, so the throw is catchable here
    // (per the class's own doc comment, doing this in a constexpr/consteval context would
    // instead be a compile error).
    EXPECT_THROW(
            (ToyEncoder{
                    std::array<std::uint8_t, 2>{0, 1},
                    std::array<std::uint8_t, 1>{99},
                    std::array<std::uint8_t, 1>{2},
            }),
            std::runtime_error);
}

TEST(test_fixed_size_encoder, qc_encoder_invalid_matrix_in_later_column_of_exponents_throws) {
    // The matrix of exponents here is 1x4 (M=1, N=4): it has only one row (row 0), so the full
    // binary parity check matrix has only expansion_factor=2 rows in total. This constructs a
    // matrix of exponents whose single nonzero entry claims row 99, which doesn't exist -- there
    // is no row 99 to shift into. If that were allowed through, encoding would compute
    // `outIdx = expansion_factor * QCrow + ... = 2*99 + ... = 198 or 199` and write there: 196+
    // positions past the end of the actual 2-row output array, i.e. an out-of-bounds write. The
    // constructor must reject this instead.
    //
    // The bad entry is deliberately placed in column 2 (of 4) of the matrix of exponents, with
    // columns 0-1 left empty, so this also confirms the check examines every column of the
    // matrix of exponents and not just the first one.
    using Enc4x2 = FixedSizeEncoderQC<std::uint8_t, /*M=*/1, /*N=*/4, /*expansion_factor=*/2, /*num_nz=*/1,
            std::uint8_t, std::uint8_t, std::uint8_t>;
    EXPECT_THROW(
            (Enc4x2{
                    std::array<std::uint8_t, 5>{0, 0, 0, 1, 1},  // colptr: only column 2 has a nonzero
                    std::array<std::uint8_t, 1>{99},             // invalid: M=1, only row 0 exists
                    std::array<std::uint8_t, 1>{0},
            }),
            std::runtime_error);
}

TEST(test_fixed_size_encoder, qc_encoder_non_power_of_two_expansion_factor) {
    // Same shape as `ToyEncoder` (a 1x1 matrix of exponents, shift=2), but expansion_factor=3, which is
    // not a power of two -- confirms the shift computation is correct beyond the power-of-two
    // codes actually used in this project's prebuilt matrices.
    using Enc3 = FixedSizeEncoderQC<std::uint8_t, /*M=*/1, /*N=*/1, /*expansion_factor=*/3, /*num_nz=*/1,
            std::uint8_t, std::uint8_t, std::uint8_t>;
    Enc3 enc{
            std::array<std::uint8_t, 2>{0, 1},
            std::array<std::uint8_t, 1>{0},
            std::array<std::uint8_t, 1>{2},  // shift = 2
    };

    // a clean cyclic shift by 2 (mod 3): in[0]->out[1], in[1]->out[2], in[2]->out[0]
    auto pos_varn = enc.get_pos_varn();
    ASSERT_EQ(pos_varn.size(), 3u);
    EXPECT_EQ(pos_varn[0], (std::vector<decltype(enc)::idx_t>{2}));
    EXPECT_EQ(pos_varn[1], (std::vector<decltype(enc)::idx_t>{0}));
    EXPECT_EQ(pos_varn[2], (std::vector<decltype(enc)::idx_t>{1}));

    std::array<std::uint8_t, 3> in{1, 0, 1};
    std::array<std::uint8_t, 3> out{};
    enc.encode_qc(std::span<std::uint8_t const, 3>{in}, std::span<std::uint8_t, 3>{out});
    EXPECT_EQ(out, (std::array<std::uint8_t, 3>{1, 1, 0}));
}
