//
// Created by Adomas Baliuka on 31.07.26.
//
// Static list of prebuilt `RateAdaptiveCode` objects
//

#ifndef LDPC4QKD_PREBUILT_CODES_HPP
#define LDPC4QKD_PREBUILT_CODES_HPP


#include <cstdint>
#include <vector>
#include <tuple>
#include <stdexcept>
#include <optional>
#include <cmath>
#include <algorithm>

#include "LDPC4QKD/fixed_size_encoder.hpp"
#include "LDPC4QKD/rate_adaptive_code.hpp"

// automatically generated LDPC codes and rate adaption arrays
#include "LDPC4QKD/autogen_ldpc_QC.hpp"
#include "LDPC4QKD/autogen/all_819k_codes.qccsc.hpp"

#include "LDPC4QKD/autogen/rate_adaption_2x4_block_4096.hpp"
#include "LDPC4QKD/autogen/rate_adaption_2x4_block_16384.hpp"
#include "LDPC4QKD/autogen/rate_adaption_2x4_block_1048576.hpp"
#include "LDPC4QKD/autogen/rate_adaption_2x6_block_6144.hpp"
#include "LDPC4QKD/autogen/rate_adaption_2x6_block_24576.hpp"
#include "LDPC4QKD/autogen/rate_adaption_2x6_block_1572864.hpp"

namespace LDPC4QKD {
    // Protograph based codes:
    constexpr inline auto encoder_2048x6144_4663d91 = helper_create_FixedSizeEncoderQC<
            AutogenLDPC_QC_2048x6144_4663d91::M, AutogenLDPC_QC_2048x6144_4663d91::expansion_factor>(
            AutogenLDPC_QC_2048x6144_4663d91::colptr, AutogenLDPC_QC_2048x6144_4663d91::row_idx,
            AutogenLDPC_QC_2048x6144_4663d91::values);

    constexpr inline auto encoder_8192x24576_71b51c1 = helper_create_FixedSizeEncoderQC<
            AutogenLDPC_QC_8192x24576_71b51c1::M, AutogenLDPC_QC_8192x24576_71b51c1::expansion_factor>(
            AutogenLDPC_QC_8192x24576_71b51c1::colptr, AutogenLDPC_QC_8192x24576_71b51c1::row_idx,
            AutogenLDPC_QC_8192x24576_71b51c1::values);

    constexpr inline auto encoder_524288x1572864_4d78a9f = helper_create_FixedSizeEncoderQC<
            AutogenLDPC_QC_524288x1572864_4d78a9f::M, AutogenLDPC_QC_524288x1572864_4d78a9f::expansion_factor>(
            AutogenLDPC_QC_524288x1572864_4d78a9f::colptr, AutogenLDPC_QC_524288x1572864_4d78a9f::row_idx,
            AutogenLDPC_QC_524288x1572864_4d78a9f::values);

    constexpr inline auto encoder_2048x4096_0c809c3 = helper_create_FixedSizeEncoderQC<
            AutogenLDPC_QC_2048x4096_0c809c3::M, AutogenLDPC_QC_2048x4096_0c809c3::expansion_factor>(
            AutogenLDPC_QC_2048x4096_0c809c3::colptr, AutogenLDPC_QC_2048x4096_0c809c3::row_idx,
            AutogenLDPC_QC_2048x4096_0c809c3::values);

    constexpr inline auto encoder_8192x16384_3fcad37 = helper_create_FixedSizeEncoderQC<
            AutogenLDPC_QC_8192x16384_3fcad37::M, AutogenLDPC_QC_8192x16384_3fcad37::expansion_factor>(
            AutogenLDPC_QC_8192x16384_3fcad37::colptr, AutogenLDPC_QC_8192x16384_3fcad37::row_idx,
            AutogenLDPC_QC_8192x16384_3fcad37::values);

    constexpr inline auto encoder_524288x1048576_9b50f98 = helper_create_FixedSizeEncoderQC<
            AutogenLDPC_QC_524288x1048576_9b50f98::M, AutogenLDPC_QC_524288x1048576_9b50f98::expansion_factor>(
            AutogenLDPC_QC_524288x1048576_9b50f98::colptr, AutogenLDPC_QC_524288x1048576_9b50f98::row_idx,
            AutogenLDPC_QC_524288x1048576_9b50f98::values);

    // degree-distribution-based codes:
    constexpr inline auto encoder_lrate_P15_block_819k = helper_create_FixedSizeEncoderQC<
            AutogenLDPC_QC_lrate_P15_block_819k::M, AutogenLDPC_QC_lrate_P15_block_819k::expansion_factor>(
            AutogenLDPC_QC_lrate_P15_block_819k::colptr, AutogenLDPC_QC_lrate_P15_block_819k::row_idx,
            AutogenLDPC_QC_lrate_P15_block_819k::values);

    constexpr inline auto encoder_lrate_P1_block_819k = helper_create_FixedSizeEncoderQC<
            AutogenLDPC_QC_lrate_P1_block_819k::M, AutogenLDPC_QC_lrate_P1_block_819k::expansion_factor>(
            AutogenLDPC_QC_lrate_P1_block_819k::colptr, AutogenLDPC_QC_lrate_P1_block_819k::row_idx,
            AutogenLDPC_QC_lrate_P1_block_819k::values);

    constexpr inline auto encoder_lrate_P25_block_819k = helper_create_FixedSizeEncoderQC<
            AutogenLDPC_QC_lrate_P25_block_819k::M, AutogenLDPC_QC_lrate_P25_block_819k::expansion_factor>(
            AutogenLDPC_QC_lrate_P25_block_819k::colptr, AutogenLDPC_QC_lrate_P25_block_819k::row_idx,
            AutogenLDPC_QC_lrate_P25_block_819k::values);

    constexpr inline auto encoder_lrate_P2_block_819k = helper_create_FixedSizeEncoderQC<
            AutogenLDPC_QC_lrate_P2_block_819k::M, AutogenLDPC_QC_lrate_P2_block_819k::expansion_factor>(
            AutogenLDPC_QC_lrate_P2_block_819k::colptr, AutogenLDPC_QC_lrate_P2_block_819k::row_idx,
            AutogenLDPC_QC_lrate_P2_block_819k::values);

    constexpr inline auto encoder_lrate_P35_block_819k = helper_create_FixedSizeEncoderQC<
            AutogenLDPC_QC_lrate_P35_block_819k::M, AutogenLDPC_QC_lrate_P35_block_819k::expansion_factor>(
            AutogenLDPC_QC_lrate_P35_block_819k::colptr, AutogenLDPC_QC_lrate_P35_block_819k::row_idx,
            AutogenLDPC_QC_lrate_P35_block_819k::values);

    constexpr inline auto encoder_lrate_P3_block_819k = helper_create_FixedSizeEncoderQC<
            AutogenLDPC_QC_lrate_P3_block_819k::M, AutogenLDPC_QC_lrate_P3_block_819k::expansion_factor>(
            AutogenLDPC_QC_lrate_P3_block_819k::colptr, AutogenLDPC_QC_lrate_P3_block_819k::row_idx,
            AutogenLDPC_QC_lrate_P3_block_819k::values);

    constexpr inline auto encoder_lrate_P45_block_819k = helper_create_FixedSizeEncoderQC<
            AutogenLDPC_QC_lrate_P45_block_819k::M, AutogenLDPC_QC_lrate_P45_block_819k::expansion_factor>(
            AutogenLDPC_QC_lrate_P45_block_819k::colptr, AutogenLDPC_QC_lrate_P45_block_819k::row_idx,
            AutogenLDPC_QC_lrate_P45_block_819k::values);

    constexpr inline auto encoder_lrate_P4_block_819k = helper_create_FixedSizeEncoderQC<
            AutogenLDPC_QC_lrate_P4_block_819k::M, AutogenLDPC_QC_lrate_P4_block_819k::expansion_factor>(
            AutogenLDPC_QC_lrate_P4_block_819k::colptr, AutogenLDPC_QC_lrate_P4_block_819k::row_idx,
            AutogenLDPC_QC_lrate_P4_block_819k::values);

    constexpr inline auto encoder_lrate_P5_block_819k = helper_create_FixedSizeEncoderQC<
            AutogenLDPC_QC_lrate_P5_block_819k::M, AutogenLDPC_QC_lrate_P5_block_819k::expansion_factor>(
            AutogenLDPC_QC_lrate_P5_block_819k::colptr, AutogenLDPC_QC_lrate_P5_block_819k::row_idx,
            AutogenLDPC_QC_lrate_P5_block_819k::values);

    constexpr inline std::tuple all_encoders_tuple{
            // Protograph-based
            encoder_2048x6144_4663d91,
            encoder_8192x24576_71b51c1,
            encoder_524288x1572864_4d78a9f,
            encoder_2048x4096_0c809c3,
            encoder_8192x16384_3fcad37,
            encoder_524288x1048576_9b50f98,
            // degree-distribution-based
            encoder_lrate_P1_block_819k,
            encoder_lrate_P15_block_819k,
            encoder_lrate_P2_block_819k,
            encoder_lrate_P25_block_819k,
            encoder_lrate_P3_block_819k,
            encoder_lrate_P35_block_819k,
            encoder_lrate_P4_block_819k,
            encoder_lrate_P45_block_819k,
            encoder_lrate_P5_block_819k
    };

    namespace HelperFixedSize {
        using Idx = std::uint32_t;
        using Bit = std::uint8_t;

        using ErrorCorrector = LDPC4QKD::RateAdaptiveCode<Idx>;

        template<typename Tout, typename Tin>
        std::vector<std::vector<Tout>> static_cast_vec_vec(std::vector<std::vector<Tin>> const &v) {
            std::vector<std::vector<Tout>> result(v.size());
            for (std::size_t i = 0; i < v.size(); ++i) {
                result[i].reserve(v[i].size());
                for (const auto &val: v[i]) {
                    result[i].push_back(val);
                }
            }
            return result;
        }

        template<typename Tout, typename Tin>
        std::vector<Tout> static_cast_vec(std::vector<Tin> const &v) {
            std::vector<Tout> result;
            result.reserve(v.size());
            for (const auto &val: v) {
                result.push_back(val);
            }
            return result;
        }

        template<typename T, std::size_t N>
        std::vector<T> arr_to_vec(std::array<T, N> a) {
            return std::vector<T>(a.begin(), a.end());
        }

        template <typename Idx=Idx>
        LDPC4QKD::RateAdaptiveCode<Idx> get_rate_adaptive_code(std::size_t id) {
            switch (id) {
                case 0: { // encoder_2048x6144_4663d91
                    auto pos_varn = static_cast_vec_vec<Idx>(
                            std::get<0>(all_encoders_tuple).get_pos_varn());
                    auto rate_adapt_rows = static_cast_vec<Idx>(
                            arr_to_vec(AutogenRateAdapt_2x6_block_6144::rows));
                    return {pos_varn, rate_adapt_rows, 0};
                }
                case 1: { // encoder_8192x24576_71b51c1
                    auto pos_varn = static_cast_vec_vec<Idx>(
                            std::get<1>(all_encoders_tuple).get_pos_varn());
                    auto rate_adapt_rows = static_cast_vec<Idx>(
                            arr_to_vec(AutogenRateAdapt_2x6_block_24576::rows));
                    return {pos_varn, rate_adapt_rows, 0};
                }
                case 2: { // encoder_524288x1572864_4d78a9f
                    auto pos_varn = static_cast_vec_vec<Idx>(
                            std::get<2>(all_encoders_tuple).get_pos_varn());
                    auto rate_adapt_rows = static_cast_vec<Idx>(
                            arr_to_vec(AutogenRateAdapt_2x6_block_1572864::rows));
                    return {pos_varn, rate_adapt_rows, 0};
                }
                case 3: { // encoder_2048x4096_0c809c3
                    auto pos_varn = static_cast_vec_vec<Idx>(
                            std::get<3>(all_encoders_tuple).get_pos_varn());
                    auto rate_adapt_rows = static_cast_vec<Idx>(
                            arr_to_vec(AutogenRateAdapt_2x4_block_4096::rows));
                    return {pos_varn, rate_adapt_rows, 0};
                }
                case 4: { // encoder_8192x16384_3fcad37
                    auto pos_varn = static_cast_vec_vec<Idx>(
                            std::get<4>(all_encoders_tuple).get_pos_varn());
                    auto rate_adapt_rows = static_cast_vec<Idx>(
                            arr_to_vec(AutogenRateAdapt_2x4_block_16384::rows));
                    return {pos_varn, rate_adapt_rows, 0};
                }
                case 5: { // encoder_524288x1048576_9b50f98
                    auto pos_varn = static_cast_vec_vec<Idx>(
                            std::get<5>(all_encoders_tuple).get_pos_varn());
                    auto rate_adapt_rows = static_cast_vec<Idx>(
                            arr_to_vec(AutogenRateAdapt_2x4_block_1048576::rows));
                    return {pos_varn, rate_adapt_rows, 0};
                }
                case 6: { // encoder 819k, lrate 0.1
                    auto pos_varn = static_cast_vec_vec<Idx>(
                            std::get<6>(all_encoders_tuple).get_pos_varn());
                    // empty: use pseudo-random (LCG-based) rows for rate adaption
                    auto rate_adapt_rows = std::vector<Idx>{};
                    return {pos_varn, rate_adapt_rows, 0};
                }
                case 7: { // encoder 819k, lrate 0.15
                    auto pos_varn = static_cast_vec_vec<Idx>(
                            std::get<7>(all_encoders_tuple).get_pos_varn());
                    // empty: use pseudo-random (LCG-based) rows for rate adaption
                    auto rate_adapt_rows = std::vector<Idx>{};
                    return {pos_varn, rate_adapt_rows, 0};
                }
                case 8: { // encoder 819k, lrate 0.2
                    auto pos_varn = static_cast_vec_vec<Idx>(
                            std::get<8>(all_encoders_tuple).get_pos_varn());
                    // empty: use pseudo-random (LCG-based) rows for rate adaption
                    auto rate_adapt_rows = std::vector<Idx>{};
                    return {pos_varn, rate_adapt_rows, 0};
                }
                case 9: { // encoder 819k, lrate 0.25
                    auto pos_varn = static_cast_vec_vec<Idx>(
                            std::get<9>(all_encoders_tuple).get_pos_varn());
                    // empty: use pseudo-random (LCG-based) rows for rate adaption
                    auto rate_adapt_rows = std::vector<Idx>{};
                    return {pos_varn, rate_adapt_rows, 0};
                }
                case 10: { // encoder 819k, lrate 0.3
                    auto pos_varn = static_cast_vec_vec<Idx>(
                            std::get<10>(all_encoders_tuple).get_pos_varn());
                    // empty: use pseudo-random (LCG-based) rows for rate adaption
                    auto rate_adapt_rows = std::vector<Idx>{};
                    return {pos_varn, rate_adapt_rows, 0};
                }
                case 11: { // encoder 819k, lrate 0.35
                    auto pos_varn = static_cast_vec_vec<Idx>(
                            std::get<11>(all_encoders_tuple).get_pos_varn());
                    // empty: use pseudo-random (LCG-based) rows for rate adaption
                    auto rate_adapt_rows = std::vector<Idx>{};
                    return {pos_varn, rate_adapt_rows, 0};
                }
                case 12: { // encoder 819k, lrate 0.4
                    auto pos_varn = static_cast_vec_vec<Idx>(
                            std::get<12>(all_encoders_tuple).get_pos_varn());
                    // empty: use pseudo-random (LCG-based) rows for rate adaption
                    auto rate_adapt_rows = std::vector<Idx>{};
                    return {pos_varn, rate_adapt_rows, 0};
                }
                case 13: { // encoder 819k, lrate 0.45
                    auto pos_varn = static_cast_vec_vec<Idx>(
                            std::get<13>(all_encoders_tuple).get_pos_varn());
                    // empty: use pseudo-random (LCG-based) rows for rate adaption
                    auto rate_adapt_rows = std::vector<Idx>{};
                    return {pos_varn, rate_adapt_rows, 0};
                }
                case 14: { // encoder 819k, lrate 0.5
                    auto pos_varn = static_cast_vec_vec<Idx>(
                            std::get<14>(all_encoders_tuple).get_pos_varn());
                    // empty: use pseudo-random (LCG-based) rows for rate adaption
                    auto rate_adapt_rows = std::vector<Idx>{};
                    return {pos_varn, rate_adapt_rows, 0};
                }
                default: {
                    throw std::runtime_error("No code available for requested ID " + std::to_string(id));
                }
            }
        }
    }

    //! Encodes `key` using the LDPC code specified by a runtime `code_id`; result is the syndrome.
    //! If `code_id` is known at compile time, prefer `encode_with_static` instead.
    //!
    //! Containers are converted to `std::span` internally: sizes are checked at runtime and must match
    //! exactly, or an exception is thrown. Containers with a compile-time-known size may instead fail to
    //! *compile*, since that size would have to match every registered code at once.
    //! Use `get_input_size`/`get_output_size` to size a `std::vector` buffer correctly, or pass `std::span`
    //! to skip the check entirely (then you must ensure the sizes are correct yourself).
    //!
    //! \tparam N internal implementation detail, need not use.
    //! \param code_id integer index into tuple of codes. Make sure both sides agree on these!
    //! \param key Contiguous container (e.g. `std::vector`, `std::array`, `std::span`) of bits (e.g. `bool` or `uint8_t`).
    //! \param result Contiguous container of bits. Used to store syndrome; must already be sized correctly for the given code!
    template<std::size_t N = 0>
    void encode_with(std::size_t code_id, auto const &key, auto &result) {
        if (code_id >= std::tuple_size_v<decltype(all_encoders_tuple)>) {
            throw std::runtime_error("Invalid code ID requested!");
        }
        if (N == code_id) {
            // if/when `all_encoders_tuple` contains non-QC matrices, this needs to change!
            // Originally, this was using `encode` instead of `encode_qc` but then it doesn't work with `vector<bool>`
            std::get<N>(all_encoders_tuple).encode_qc(key, result);
            return;
        }

        if constexpr (N + 1 < std::tuple_size_v<decltype(all_encoders_tuple)>) {
            return encode_with<N + 1>(code_id, key, result);
        }
    }

    //! Same as `encode_with`, but for a `code_id` known at compile time.
    //! \tparam code_id integer index into tuple of codes. Make sure both sides agree on these!
    //! \param key Contiguous container (e.g. `std::vector`, `std::array`, `std::span`) of bits (e.g. `bool` or `uint8_t`).
    //! \param result Contiguous container of bits. Used to store syndrome; must already be sized correctly for the given code!
    template<std::size_t code_id>
    void encode_with_static(auto const &key, auto &result) {
        std::get<code_id>(all_encoders_tuple).encode(key, result);
    }

    //! Get input size of code with given ID.
    //!
    //! \tparam N internal implementation detail, need not use.
    //! \param code_id integer index into tuple of codes. Make sure both sides agree on these!
    template<std::size_t N = 0>
    constexpr std::size_t get_input_size(std::size_t code_id) {
        if (code_id >= std::tuple_size_v<decltype(all_encoders_tuple)>) {
            throw std::runtime_error("Invalid code ID requested!");
        }
        if (N == code_id) {
            return std::get<N>(all_encoders_tuple).get_input_size();
        }

        if constexpr (N + 1 < std::tuple_size_v<decltype(all_encoders_tuple)>) {
            return get_input_size<N + 1>(code_id);
        }
        return 0; // this should never happen
    }

    //! Get output size (number rows of parity check matrix) of code with given ID.
    //!
    //! \tparam N internal implementation detail, need not use.
    //! \param code_id integer index into tuple of codes. Make sure both sides agree on these!
    template<std::size_t N = 0>
    constexpr std::size_t get_output_size(std::size_t code_id) {
        if (code_id >= std::tuple_size_v<decltype(all_encoders_tuple)>) {
            throw std::runtime_error("Invalid code ID requested!");
        }
        if (N == code_id) {
            return std::get<N>(all_encoders_tuple).get_output_size();
        }

        if constexpr (N + 1 < std::tuple_size_v<decltype(all_encoders_tuple)>) {
            return get_output_size<N + 1>(code_id);
        }
        return 0; // this should never happen
    }

    //! Shannon binary entropy function (in bits).
    //! \param p must be in [0, 1].
    inline double binary_entropy(double p) {
        if (p < 0 || p > 1) {
            throw std::domain_error("binary_entropy: p must be between 0 and 1");
        }
        if (p == 0 || p == 1) {
            return 0;
        }
        return -p * std::log2(p) - (1 - p) * std::log2(1 - p);
    }

    //! Result of `select_suitable_code`: a prebuilt code ID and the recommended syndrome size (in bits) for a
    //! single block of that code.
    struct SuitableCodeChoice {
        std::size_t code_id;
        std::size_t syndrome_bits_per_block;
    };

    //! Select a suitable prebuilt LDPC code and a recommended syndrome size for a given estimated channel
    //! parameter.
    //!
    //! \param ch_param_estimate estimated bit-flip probability of a binary symmetric channel.
    //! \param input_block_size size (in bits) of the caller's input block.
    //! \return chosen code ID and suggested syndrome size, or `std::nullopt` if no code available.
    //! The selected code may have a **lower** `input_block_size` than requested but not higher.
    inline std::optional<SuitableCodeChoice> select_suitable_code(
            double ch_param_estimate,
            [[maybe_unused]] std::size_t input_block_size) {
        if (ch_param_estimate <= 0) {
            throw std::invalid_argument("ch_param_estimate estimate must be > 0");
        }

        std::size_t code_id;
        double target_lrate;
        if (ch_param_estimate < 0.01) {
            code_id = 1;
            target_lrate = 1. / 6.;
        } else if (ch_param_estimate < 0.03) {
            code_id = 1;
            target_lrate = 3. * binary_entropy(ch_param_estimate);
        } else if (ch_param_estimate < 0.049) {
            code_id = 1;
            target_lrate = 1.8 * binary_entropy(ch_param_estimate);
        } else if (ch_param_estimate < 0.07) {
            code_id = 4;
            target_lrate = 1.7 * binary_entropy(ch_param_estimate);
        } else if (ch_param_estimate < 0.092) {
            code_id = 4;
            target_lrate = 1.3 * binary_entropy(ch_param_estimate);
        } else {
            return std::nullopt;
        }

        const auto code = HelperFixedSize::get_rate_adaptive_code(code_id);
        const auto syndrome_bits_per_block = std::min<std::size_t>(
                code.get_n_rows_mother_matrix(),
                static_cast<std::size_t>(std::floor(static_cast<double>(code.getNCols()) * target_lrate)));

        return SuitableCodeChoice{code_id, syndrome_bits_per_block};
    }

}

#endif //LDPC4QKD_PREBUILT_CODES_HPP
