//
// Created by alice on 23.04.21.
// Note: this file is for applications (e.g. small executable required) where only the encoder (NO decoder) is desired.
// The decoder can also encode and has a much cleaner interface.
// Basically, this file just implements matrix-vector multiplication in compressed sparse column (CSC) storage format.
//

#ifndef LDPC4QKD_ENCODER_HPP
#define LDPC4QKD_ENCODER_HPP

// Standard library
#include <cstdint>
#include <array>
#include <algorithm>


#ifdef DEBUG_MESSAGES_ENABLED

#include <iostream>

#define DEBUG_MESSAGE(msg) do {std::cerr << msg << std::endl;} while (false)

#else

#define DEBUG_MESSAGE(msg)

#endif /* ifdef DEBUG_MESSAGES_ENABLED */

// automatically generated code containing the LDPC matrix and
// the line indices to be combined for rate adaption.
#include "autogen_ldpc_matrix_csc.hpp"
#include "autogen_rate_adaption.hpp"


namespace LDPC4QKD {
    template<typename Bit>
    constexpr bool xor_as_bools(Bit lhs, Bit rhs) {
        return (static_cast<bool>(lhs) != static_cast<bool>(rhs));
    }


    template<typename Bit>
    constexpr void encode(const std::array<Bit, AutogenLDPC::N> &in, std::array<Bit, AutogenLDPC::M> &out) {
        using namespace AutogenLDPC;
        static_assert(AutogenLDPC::N >= AutogenLDPC::M, "The syndrome should be shorter than the input bitstring.");

        for (std::size_t col = 0; col < in.size(); col++) {
            for (std::size_t j = colptr[col]; j < colptr[col + 1]; j++) {
                out[row_idx[j]] = xor_as_bools(out[row_idx[j]], in[col]);
            }
        }
    }


    /// TODO reorder rows?
    // This is just for testing and seeing the difference in performance
    // The rate adaptive matrix can be defined in Julia as
    //```julia
    //function simple_rate_adaption(
    //    H::AbstractArray{Int8, 2},
    //    ra_lines::AbstractArray{T, 2} where T <: Integer,  # like C14
    //    n_line_combs::Integer
    //    )
    //    kpos = setdiff(1:size(H, 1), ra_lines[1:n_line_combs, :])
    //    Hra = [H[kpos, :] ; rem.(H[ra_lines[1:n_line_combs, 1], :] + H[ra_lines[1:n_line_combs, 2], :], 2)]
    //
    //    return Hra
    //end
    //```
    template<typename Bit, std::size_t reduced_size>
    constexpr void rate_adapt(
            const std::array<Bit, AutogenLDPC::M> &syndrome, std::array<Bit, reduced_size> &reduced_syndrome) {
        // Depending on the rate adaption requested, only part of the `rows` array will be used.
        using AutogenRateAdapt::rows;

        // have to modify values to mark which bits are used during rate adaption
        // Also need -1 as a value, which has special purpose below.
        std::array<signed short, AutogenLDPC::M> syndrome_copy;  // uninitialized because filled via copy right below
        std::copy(syndrome.begin(), syndrome.end(), syndrome_copy.begin());

        static_assert(reduced_size < AutogenLDPC::M,
                "Requested rate adapted syndrome size must be less than the original syndrome size.");

        constexpr std::size_t n_row_combinations = AutogenLDPC::M - reduced_size;
        static_assert(rows.size() / 2 >= n_row_combinations,
                "The specified rate adaption does not support such a high amount of line combinations.");

        std::size_t start_of_ra_part = reduced_syndrome.size() - n_row_combinations;

        // put results of combined lines at the back of output.
        for (std::size_t i{}; i < n_row_combinations; ++i) {
            reduced_syndrome[start_of_ra_part + i] = xor_as_bools(syndrome_copy[rows[2 * i]],
                                                                  syndrome_copy[rows[2 * i + 1]]);
            syndrome_copy[rows[2 * i]] = -1;  // -1 marks that the value has been used.
            syndrome_copy[rows[2 * i + 1]] = -1;
        }

        std::size_t j{};
        // put the remaining bits that were not rate adapted at the front of output.
        for (std::size_t i{}; i < start_of_ra_part; ++i) {
            while (syndrome_copy[j] == -1) {
                j++;
            }
            reduced_syndrome[i] = syndrome_copy[j];
            j++;
        }
    }


    template<typename Bit>
    [[deprecated("TODO this function is probably slow and should not be used. Update to be like 'safe' version?")]]
    constexpr void rate_adapt_unsafe(
            const std::array<Bit, AutogenLDPC::M> &syndrome,
            Bit *rate_adapted_syndrome,
            const std::size_t size_rate_adapted_syndrome) {
        using AutogenRateAdapt::rows;
        if(size_rate_adapted_syndrome >= syndrome.size()) {
            DEBUG_MESSAGE("Requested rate adapted syndrome size must be less than the original syndrome size.");
            // TODO consider removing exception?
            throw std::runtime_error("Requested rate adapted syndrome size must be less than the original syndrome size.");
        }
        std::size_t n_row_combs = syndrome.size() - size_rate_adapted_syndrome;
        if(rows.size() / 2 < n_row_combs) {
            DEBUG_MESSAGE("The specified rate adaption does not support such a high amount of line combinations");
            // TODO consider removing exception?
            throw std::runtime_error("The specified rate adaption does not support such a high amount of line combinations");
        }

        // Depending on the rate adaption requested, only part of the `rows` array will be used.
        // This corresponds to `used_rows.end()`,
        // where `used_rows` is the portion of `rows` that is relevant for the specified rate (given by `reduced_size`)
        auto rows_end_indx = rows.begin() + 2 * n_row_combs;

        std::size_t reduced_syndrome_idx{};
        for (std::size_t i{}; i < syndrome.size(); ++i) {
            // if `i` is not the list of line indices to be combined,
            // add the corresponding syndrome bit to the `reduced_syndrome`
            if (std::find(rows.begin(), rows_end_indx, i) == rows_end_indx) {
                rate_adapted_syndrome[reduced_syndrome_idx] = syndrome[i];
                reduced_syndrome_idx++;
            }
        }

        // combine the lines as specified by the indices in the array `AutogenRateAdapt::rows`.
        for (std::size_t i{}; i < 2 * n_row_combs; i += 2) {
            rate_adapted_syndrome[reduced_syndrome_idx] = xor_as_bools(syndrome[rows[i]], syndrome[rows[i + 1]]);
            reduced_syndrome_idx++;
        }
    }


    /// TODO in-place computation may require moving bits around a lot. Think about how to achieve it efficiently.
//    /// TODO allow runtime determination of reduced_size!!!
//    ///     This is just for testing and seeing the difference in performance
//    template<typename Bit, std::size_t reduced_size>
//    constexpr void rate_adapt_inplace(const std::array<Bit, AutogenLDPC::M> &syndrome) {
//
//    }
//
//    /// TODO allow runtime determination of reduced_size!!!
//    ///     This is just for testing and seeing the difference in performance
//    template<typename Bit>
//    constexpr void rate_adapt_inplace(const std::array<Bit, AutogenLDPC::M> &syndrome, const std::size_t reduced_size) {
//
//    }

}  // namespace RALDPC

#endif //LDPC4QKD_ENCODER_HPP
