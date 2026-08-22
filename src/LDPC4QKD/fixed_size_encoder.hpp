//
// Created by Adomas Baliuka on 03.12.24.
//
// Representation of fixed-size (dimensions known at compile time) LDPC encoder (syndrome computation).
// Supports quasi-cyclic (QC) format as input.
//

#ifndef QKD_POSTPROC_BOB_ENCODER_ADVANCED_HPP
#define QKD_POSTPROC_BOB_ENCODER_ADVANCED_HPP

#include <cstdint>
#include <array>
#include <span>
#include <vector>
#include <concepts>
#include <sstream>
#include <stdexcept>

namespace LDPC4QKD {
    template<std::size_t N>
    constexpr auto bits_needed() {
        auto n = N;
        std::size_t number_of_bits{0ul};
        while (n > 0) {
            n >>= 1;
            number_of_bits++;
        }
        return number_of_bits;
    }

    template<std::size_t N>
    using smallest_type = std::conditional_t<bits_needed<N>() <= 8 * sizeof(std::uint8_t), std::uint8_t,
            std::conditional_t<bits_needed<N>() <= 8 * sizeof(std::uint16_t), std::uint16_t,
                    std::conditional_t<bits_needed<N>() <= 8 * sizeof(std::uint32_t), std::uint32_t, std::uint64_t>>>;


    template<typename Bit1, typename Bit2>
    constexpr bool xor_as_bools(Bit1 lhs, Bit2 rhs) {
        return (static_cast<bool>(lhs) != static_cast<bool>(rhs));
    }

    struct FixedSizeInputOutput {
        [[nodiscard]] virtual constexpr std::size_t get_input_size() const = 0;

        [[nodiscard]] virtual constexpr std::size_t get_output_size() const = 0;

        virtual constexpr ~FixedSizeInputOutput() = default;
    };

    template<typename idx_t>
    struct ComputablePosVar : public FixedSizeInputOutput {
        [[nodiscard]] virtual constexpr std::vector<std::vector<idx_t>> get_pos_varn() const = 0;

        constexpr ~ComputablePosVar() override = default;
    };


    template<
            typename bit_type, /// e.g. bool, or std::uint8_t
            std::size_t output_size, std::size_t input_size>
    struct FixedSizeEncoder : public ComputablePosVar<smallest_type<input_size>> {
        static constexpr std::size_t outputSize = output_size;
        static constexpr std::size_t inputSize = input_size;

        [[nodiscard]] constexpr std::size_t get_input_size() const override {
            return inputSize;
        }

        [[nodiscard]] constexpr std::size_t get_output_size() const override {
            return output_size;
        }

        /// performant but no size check (user has to provide valid std::span)
        /// key shall not be changed!
        /// also: inputs/outputs have to be contiguous in memory.
        virtual constexpr void encode_span(
                std::span<bit_type const, input_size> key,
                std::span<bit_type, output_size> syndrome) const = 0;

        /// general syndrome computation with runtime size check (a runtime cost usually not worth worrying about).
        /// Not constexpr: using `std::stringstream` in constexpr context only possible in C++23
        //  Compile-time checks use `encode_span` instead.
        void encode(auto const &key, auto &syndrome) const {
            if (key.size() == input_size && syndrome.size() == output_size) {
                encode_span(
                        std::span<bit_type const, input_size>{key},
                        std::span<bit_type, output_size>{syndrome});
            } else {
                std::stringstream s;
                s << "LDPC encoder: incorrect sizes of intput / output arrays\n"
                  << "RECEIVED: key.size() = " << key.size() << ". " << "syndrome.size() = " << syndrome.size() << ".\n"
                  << "EXPECTED: key.size() = " << input_size << ". " << "syndrome.size() = " << output_size << ".\n";
                throw std::out_of_range(s.str());
            }
        }

        /// Note: vector<bool> cannot be supported in the interface specified here.
        /// This is due to the template specialization for vector<bool>,
        /// which does not allow direct conversion to span, as in `std::span<bool, N>{vec}`.
        void encode(std::vector<bool> const &key, std::vector<bool> &syndrome) const = delete;

        // Have to write explicitly for some gcc versions. See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=93413
        constexpr ~FixedSizeEncoder() override = default;
    };

    template<
            typename bit_type, /// e.g. bool, or std::uint8_t
            std::size_t M, std::size_t N, std::size_t expansion_factor, std::size_t num_nz,
            std::integral coptr_uintx_t, std::integral row_idx_uintx_t, std::integral values_uintx_t>
    struct FixedSizeEncoderQC : public FixedSizeEncoder<bit_type, M * expansion_factor, N * expansion_factor> {

        constexpr FixedSizeEncoderQC(std::array<coptr_uintx_t, N + 1> colptr,
                                     std::array<row_idx_uintx_t, num_nz> row_idx,
                                     std::array<values_uintx_t, num_nz> values) :
                colptr(colptr), row_idx(row_idx), values(values) {
            if (!matrix_consistent_with_input_size()) {
                // If constructed at compile time, this throw surfaces as "expression is not a constant
                // expression" rather than showing the message.
                // If that happens, assume there is an issue with the given inputs.
                throw std::runtime_error("Inputs would make encoder that performs out-of-memory access!");
            }
        }

        constexpr void encode_span(
                std::span<bit_type const, N * expansion_factor> key,
                std::span<bit_type, M * expansion_factor> syndrome) const override {
            encode_qc(key, syndrome);
        }

        // Have to write explicitly for some gcc versions. See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=93413
        constexpr ~FixedSizeEncoderQC() override = default;

        using idx_t = smallest_type<FixedSizeEncoderQC::inputSize>;

        [[nodiscard]] constexpr std::vector<std::vector<idx_t>> get_pos_varn() const override {
            const auto n_cols = FixedSizeEncoderQC::inputSize;
            const auto n_rows = FixedSizeEncoderQC::outputSize;

            std::vector<std::vector<idx_t>> pos_varn{};
            pos_varn.assign(n_rows, std::vector<idx_t>{});

            for (idx_t col = 0; col < n_cols; col++) {
                auto QCcol = col / expansion_factor;  // column index into matrix of exponents
                auto col_in_block = (col % expansion_factor);  // column's position within its own QC block
                for (auto j = colptr[QCcol]; j < colptr[QCcol + 1]; j++) {
                    auto shiftVal = values[j];
                    auto QCrow = row_idx[j];  // row index into matrix of exponents
                    // `outIdx` is the row (into the full matrix) at which this sub-block places a `1`:
                    // the sub-block's base row `expansion_factor * QCrow`, plus the shifted position of
                    // `col` within the block.
                    // `col_in_block` and `shift_in_block` are both already in [0, expansion_factor),
                    // so subtracting and re-adding `expansion_factor` before the final `% expansion_factor`
                    // keeps the result correct (and non-negative) for any expansion_factor.
                    auto shift_in_block = (shiftVal % expansion_factor);
                    auto outIdx = (expansion_factor * QCrow) +
                        ((col_in_block + expansion_factor - shift_in_block) % expansion_factor);
                    pos_varn[outIdx].push_back(static_cast<idx_t>(col));
                }
            }

            return pos_varn;
        }

        /// Avoiding runtime length-check from the types.
        /// TODO this template overload does not match things like `std::array`, although it would be nice!
        constexpr void encode_qc(
                std::span<bit_type const, N * expansion_factor> in,
                std::span<bit_type, M * expansion_factor> out) const {

            static_assert(N >= M, "The syndrome should be shorter than the input bitstring.");

            for (std::size_t col = 0; col < in.size(); col++) {
                auto QCcol = col / expansion_factor;  // column index into matrix of exponents
                auto col_in_block = (col % expansion_factor);  // column's position within its own QC block
                for (std::size_t j = colptr[QCcol]; j < colptr[QCcol + 1]; j++) {
                    auto shiftVal = values[j];
                    auto QCrow = row_idx[j];  // row index into matrix of exponents
                    // outIdx formula: see the comment in get_pos_varn() above.
                    auto shift_in_block = (shiftVal % expansion_factor);
                    auto outIdx = (expansion_factor * QCrow) +
                        ((col_in_block + expansion_factor - shift_in_block) % expansion_factor);

                    out[outIdx] = xor_as_bools(out[outIdx], in[col]);
                }
            }
        }

        /// General overload, which does not assume spans, but does a **runtime length check!**.
        /// Note: deliberately NOT constexpr, for the same reason as `FixedSizeEncoder::encode` above.
        void encode_qc(auto const &in, auto &out) const {
            if (std::size(in) != N * expansion_factor || std::size(out) != M * expansion_factor) {
                std::stringstream s;
                s << "LDPC encoder: incorrect sizes of intput / output arrays\n"
                  << "RECEIVED: key.size() = " << in.size() << ". " << "syndrome.size() = " << out.size() << ".\n"
                  << "EXPECTED: key.size() = " << N * expansion_factor << ". "
                  << "syndrome.size() = " << M * expansion_factor << ".\n";
                throw std::out_of_range(s.str());
            }

            for (std::size_t col = 0; col < in.size(); col++) {
                auto QCcol = col / expansion_factor;  // column index into matrix of exponents
                auto col_in_block = (col % expansion_factor);  // column's position within its own QC block
                for (std::size_t j = colptr[QCcol]; j < colptr[QCcol + 1]; j++) {
                    auto shiftVal = values[j];
                    auto QCrow = row_idx[j];  // row index into matrix of exponents
                    // outIdx formula: see the comment in get_pos_varn() above.
                    auto shift_in_block = (shiftVal % expansion_factor);
                    auto outIdx = (expansion_factor * QCrow) +
                        ((col_in_block + expansion_factor - shift_in_block) % expansion_factor);

                    out[outIdx] = xor_as_bools(out[outIdx], in[col]);
                }
            }
        }

    private:
        /// Checks that this QC-encoder never accesses input/output arrays out of bounds.
        /// Runs at compile time for every prebuilt code.
        /// If false, the object is invalid and must not be used.
        /// Two conditions are checked:
        ///  1. `colptr` is non-decreasing and bounded by `num_nz` (it's used as the [begin, end) range into
        ///     `row_idx`/`values`).
        ///  2. every `row_idx[j] < M` -- since `outIdx = expansion_factor*QCrow + r` with `r` always in
        ///     `[0, expansion_factor)`, this is equivalent to `outIdx < M*expansion_factor`, so checking each
        ///     row once suffices without expanding to the full column range
        [[nodiscard]] constexpr bool matrix_consistent_with_input_size() const {
            for (std::size_t QCcol = 0; QCcol < N; QCcol++) {
                if (colptr[QCcol + 1] < colptr[QCcol] || colptr[QCcol + 1] > num_nz) {
                    return false;
                }
            }

            for (std::size_t QCcol = 0; QCcol < N; QCcol++) {
                for (std::size_t j = colptr[QCcol]; j < colptr[QCcol + 1]; j++) {
                    auto QCrow = row_idx[j];  // row index into matrix of exponents
                    if (QCrow >= M) {
                        return false;
                    }
                }
            }
            return true;
        }

        std::array<coptr_uintx_t, N + 1> colptr;
        std::array<row_idx_uintx_t, num_nz> row_idx;
        std::array<values_uintx_t, num_nz> values;
    };


    /// Deduces `bit_type`/`N`/`num_nz` from the array arguments,
    /// so callers don't have to spell out every template parameter of `FixedSizeEncoderQC` manually.
    template<std::size_t M, std::size_t expansion_factor>
    consteval auto helper_create_FixedSizeEncoderQC(auto colptr, auto row_idx, auto values) {
        using bit_type = std::uint8_t;
        constexpr auto num_nz = values.size();
        constexpr auto N = colptr.size() - 1;
        using FixedSizeEncoderQC_inst = FixedSizeEncoderQC<bit_type, M, N, expansion_factor, num_nz,
                typename decltype(colptr)::value_type,
                typename decltype(row_idx)::value_type,
                typename decltype(values)::value_type>;
        return FixedSizeEncoderQC_inst{colptr, row_idx, values};
    }

}

#endif //QKD_POSTPROC_BOB_ENCODER_ADVANCED_HPP
