/*!
 * This header file contains a belief propagation (BP) decoder for low density parity check (LDPC) codes.
 * It also supports rate adaption (reducing the number of LDPC matrix rows).
 */


#ifndef LDPC4QKD_LDPC_MATRIX_HPP
#define LDPC4QKD_LDPC_MATRIX_HPP

#include <cstdint>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <concepts>
#include <utility>
#include <limits>

#include "rate_adaption_random.hpp"

#ifdef LDPC4QKD_DEBUG_MESSAGES_ENABLED

#include <iostream>

#define LDPC4QKD_DEBUG_MESSAGE(msg) do {std::cerr << msg << std::endl;} while (false)

#else

#define LDPC4QKD_DEBUG_MESSAGE(msg)

#endif /* ifdef LDPC4QKD_DEBUG_MESSAGES_ENABLED */


namespace LDPC4QKD {

    //! Selects which belief-propagation variant `decode_at_current_rate`/`decode_infer_rate` runs.
    enum class Decoder {
        Flooding,  //!< Original flooding-schedule BP (see `decode_flooding`). Simplest decoder.
        Layered,   //!< Layered/serial-schedule BP (see `decode_layered`): default. Typically converges in about
                   //!< half the iterations of flooding, for equal or better FER.
        Improved,  //!< Layered BP with damping and a bit-flip rescue stage (see `decode_improved`).
    };

    inline double tanh_half(double x) {
        auto exp_x = ::exp(x);
        return (exp_x - 1) / (exp_x + 1);
        // Note: this seems to be faster than `::tanh(0.5 * x)`.
    }

    /*!
     * compute log-likelihood-ratios for a given keye and channel parameter
     * @tparam Bit e.g. std::uint8_t or bool
     * @param bitstring string of bits
     * @param bsc_channel_parameter bit-flip probability of binary symmetric channel (BSC)
     * @return log-likelyhoods corresponding to input bitstring
     */
    template<typename Bit>
    std::vector<double> llrs_bsc(const std::vector<Bit> &bitstring, const double bsc_channel_parameter) {
        double vlog = ::log((1 - bsc_channel_parameter) / bsc_channel_parameter);
        std::vector<double> llrs(bitstring.size());
        for (std::size_t i{}; i < llrs.size(); ++i) {
            llrs[i] = vlog * (1 - 2 * bitstring[i]); // log likelihood ratios
        }
        return llrs;
    }

    /*!
     * Belief propagation (BP) decoder for binary low density parity check (LDPC) codes.
     * Supports rate adaption (reducing the number of LDPC matrix rows).
     * Intended for distributed source coding (a.k.a. Slepian-Wolf coding).
     * LDPC code is stored in sparse column storage (CSC) format.
     *
     *
     * @tparam idx_t unsigned integer type fitting number of columns N (thus also number of rows M)
     */
    template<std::unsigned_integral idx_t=std::uint16_t>
    class RateAdaptiveCode {
    public:
        // ------------------------------------------------------------------------------------------------ type aliases
        using MatrixIndex = idx_t;

        // ------------------------------------------------------------------------------------------------ constructors
        /*!
         * Constructor for using the code without rate adaption.
         * The parity check matrix matrix is stored using Compressed Sparse Column (CSC) format.
         *
         * @tparam colptr_t unsigned integer type that fits ("number of non-zero matrix entries" + 1)
         * @param colptr column pointer array for specifying mother parity check matrix.
         * @param rowIdx row index array for specifying mother parity check matrix.
         */
        template<typename colptr_t>
        RateAdaptiveCode(const std::vector<colptr_t> &colptr, const std::vector<idx_t> &rowIdx)
                : n_mother_rows(*std::max_element(rowIdx.begin(), rowIdx.end()) + 1u),
                  n_cols(colptr.size() - 1),
                  mother_pos_varn(compute_mother_pos_varn(colptr, rowIdx)), // computed here and henceforth `const`!
                  rows_to_combine({}),
                  auto_rate_adaption(RateAdaptLCG::get_LCG_with_period(n_mother_rows, 0)) {
            constexpr idx_t n_line_combs = 0;
            recompute_pos_vn_cn(n_line_combs);
        }

        /*!
         * Constructor for using the code with rate adaption.
         * The mother parity check matrix is stored using Compressed Sparse Column (CSC) format.
         * The rate adaption is stored as an array of matrix row indices, which are combined for rate adaption.
         *
         * note: invalid `rows_to_combine_rate_adapt`, for example non-zero based, may lead to a segmentation fault.
         *
         * note: repeated node indices (variable node eliminations) after rate adaption are detected and
         *      removed automatically during `recompute_pos_vn_cn`.
         *
         * @tparam colptr_t unsigned integer type that fits ("number of non-zero matrix entries" + 1)
         * @param colptr column pointer array for specifying mother parity check matrix.
         * @param rowIdx row index array for specifying mother parity check matrix.
         * @param rows_to_combine_rate_adapt array of mother-matrix line indices to be combined for rate adaption
         * @param initial_row_combs number of line indices to combine initially
         */
        template<typename colptr_t>
        RateAdaptiveCode(std::vector<colptr_t> colptr,
                         std::vector<idx_t> rowIdx,
                         std::vector<idx_t> rows_to_combine_rate_adapt,
                         idx_t initial_row_combs = 0)
                : n_mother_rows(*std::max_element(rowIdx.begin(), rowIdx.end()) + 1u),
                  n_cols(colptr.size() - 1),
                  mother_pos_varn(compute_mother_pos_varn(colptr, rowIdx)), // computed here and henceforth `const`!
                  rows_to_combine(std::move(rows_to_combine_rate_adapt)),
                  auto_rate_adaption(RateAdaptLCG::get_LCG_with_period(n_mother_rows, 0)) {
            if (rows_to_combine.size() % 2 != 0) {
                throw std::domain_error("The number of rows to combine for rate adaption "
                                        "(size of argument array) is an odd number (expected even).");
            }

            if (initial_row_combs > get_max_ra_steps()) {
                throw std::domain_error("The number of desired initial row combinations for rate adaption "
                                        "is larger than the given array of lines to combine.");
            }

            // compute current `pos_varn` and `pos_checkn` from `mother_pos_varn`
            recompute_pos_vn_cn(initial_row_combs);
        }

        /*!
         * Constructor for using the code with rate adaption.
         * The mother parity check matrix is stored in `mother_pos_varn`
         * The rate adaption is stored as an array of matrix row indices, which are combined for rate adaption.
         *
         * @param mother_pos_varn input check nodes to each variable node of the mother matrix
         * @param rows_to_combine_rate_adapt array of mother-matrix line indices to be combined for rate adaption
         * @param initial_row_combs number of line indices to combine initially
         */
        RateAdaptiveCode(std::vector<std::vector<idx_t>> mother_pos_varn,
                         std::vector<idx_t> rows_to_combine_rate_adapt,
                         idx_t initial_row_combs = 0)
                : n_mother_rows(mother_pos_varn.size()),
                  n_cols(compute_n_cols(mother_pos_varn)),
                  mother_pos_varn(std::move(mother_pos_varn)), // computed here and henceforth `const`!
                  rows_to_combine(std::move(rows_to_combine_rate_adapt)),
                  auto_rate_adaption(RateAdaptLCG::get_LCG_with_period(n_mother_rows, 0)),
                  n_ra_rows(n_mother_rows - initial_row_combs) {
            if (rows_to_combine.size() % 2 != 0) {
                throw std::domain_error("The number of rows to combine for rate adaption "
                                        "(size of argument array) is an odd number (expected even).");
            }

            if (initial_row_combs > get_max_ra_steps()) {
                throw std::domain_error("The number of desired initial row combinations for rate adaption "
                                        "is larger than the given array of lines to combine.");
            }

            // compute current `pos_varn` and `pos_checkn` from `mother_pos_varn`
            recompute_pos_vn_cn(initial_row_combs);
        }

        // ---------------------------------------------------------------------------------------------- public methods


        /*!
         *  Encode (i.e., compute syndrome) using mother matrix
         * @tparam BitL e.g. std::uint8_t or bool.
         * @tparam BitRBitR allowed to be signed, to enable the "mark combined by -1" trick in `encode_with_ra`.
         * @param in input bitvector
         * @param out output bitvector
         */
        template<typename BitL=bool, typename BitR=bool>
        constexpr void encode_no_ra(const std::vector<BitL> &in, std::vector<BitR> &out) const {
            if (in.size() != n_cols) {
                throw std::domain_error("Encoder (encode_no_ra) received invalid input length.");
            }
            out.assign(n_mother_rows, 0);

            for (std::size_t i{}; i < mother_pos_varn.size(); ++i) {
                for (auto &var_node: mother_pos_varn[i]) {
                    out[i] = xor_as_bools(out[i], in[var_node]);
                }
            }
        }

        /*!
         * Compute syndrome using given rate adaption. Does not change internal rate adaption state!
         * @tparam Bit e.g. std::uint8_t or bool.
         * @param in input array
         * @param out Vector to store syndrome. Will be resized to `output_syndrome_length`
         * @param output_syndrome_length Desired length of syndrome (exception is thrown if not satisfiable)
         */
        template<std::unsigned_integral Bit>
        void encode_with_ra(
                const std::vector<Bit> &in, std::vector<Bit> &out, std::size_t output_syndrome_length) const {
            if (in.size() != n_cols) {
                throw std::domain_error("Encoder (encode_with_ra) received invalid input length.");
            }
            if (output_syndrome_length > n_mother_rows) {
                throw std::domain_error("Requested syndrome is larger than the number of rows of the mother matrix.");
            }
            if (output_syndrome_length < n_mother_rows - get_max_ra_steps()) {
                throw std::domain_error("Requested syndrome is smaller than supported by the specified rate adaption.");
            }

            // HAVE to use !!!SIGNED!!! `int8_t`!!! Value `-1` is used below to mark bits included into final output!
            std::vector<std::int8_t> non_ra_encoding;
            encode_no_ra(in, non_ra_encoding);

            // now use the non-rate adapted syndrome to compute the rate adapted syndrome
            const std::size_t n_line_combinations = n_mother_rows - output_syndrome_length;
            out.assign(output_syndrome_length, 0);

            // If `rows_to_combine` is empty, pairs are auto-generated (LCG-based) and placed at the FRONT of the
            // output; otherwise the explicitly given pairs are used, placed at the BACK.
            const bool auto_generated = rows_to_combine.empty();
            RateAdaptLCG::LCG local_lcg = auto_rate_adaption;

            const std::size_t start_of_ra_part = auto_generated
                    ? 0 : (output_syndrome_length - n_line_combinations);

            for (std::size_t i{}; i < n_line_combinations; ++i) {
                idx_t idx1, idx2;
                if (auto_generated) {
                    idx1 = static_cast<idx_t>(local_lcg.next());
                    idx2 = static_cast<idx_t>(local_lcg.next());
                } else {
                    idx1 = rows_to_combine[2 * i];
                    idx2 = rows_to_combine[2 * i + 1];
                }

                if (auto_generated) {
                    // LCG-selected pairs are not guaranteed variable-disjoint, so the two mother-code
                    // syndrome bits alone don't determine this check's true value.
                    // To compute the rate-adapted syndrome, we must recompute its parity
                    // directly from `in`, over the union of both rows' variable lists.
                    // This mirrors recompute_pos_vn_cn exactly, and is slower than the explicit version, when `auto_generated == false`.
                    // Also slower than `encode_at_current_rate`, however, still faster than calling both `set_rate` and `encode_at_current_rate`.
                    std::vector<idx_t> union_vars{mother_pos_varn[idx1]};
                    union_vars.insert(union_vars.end(), mother_pos_varn[idx2].begin(), mother_pos_varn[idx2].end());
                    std::sort(union_vars.begin(), union_vars.end());
                    union_vars.erase(std::unique(union_vars.begin(), union_vars.end()), union_vars.end());

                    Bit combined = 0;
                    for (auto v : union_vars) {
                        combined = xor_as_bools(combined, in[v]);
                    }
                    out[start_of_ra_part + i] = combined;
                } else {
                    // Explicit pairs are pre-verified variable-disjoint, so XOR of the two already-computed
                    // mother-code syndrome bits is exact and cheaper than recomputing from `in`.
                    out[start_of_ra_part + i] = xor_as_bools(non_ra_encoding[idx1], non_ra_encoding[idx2]);
                }
                non_ra_encoding[idx1] = -1;  // -1 marks that the value has been used.
                non_ra_encoding[idx2] = -1;
            }

            std::size_t j{};
            // put the remaining bits that were not rate adapted into the rest of output
            // (after the combined bits for auto-generated rate adaption; before them otherwise).
            const std::size_t leftover_begin = auto_generated ? n_line_combinations : 0;
            const std::size_t leftover_end = auto_generated ? output_syndrome_length : start_of_ra_part;
            for (std::size_t i = leftover_begin; i < leftover_end; ++i) {
                while (non_ra_encoding[j] == -1) {
                    j++;
                }
                out[i] = non_ra_encoding[j];
                j++;
            }
        }

        /// decoder infers rate from the length of the syndrome and changes the internal decoder state to match this rate.
        /// Note: since this function modifies the code (by performing rate adaption), it is NOT CONST.
        /// this change may be somewhat computationally expensive
        /// `Bit` should be e.g. std::uint8_t or bool.
        template<std::unsigned_integral Bit>
        bool decode_infer_rate(const std::vector<double> &llrs,
                               const std::vector<Bit> &syndrome,
                               std::vector<Bit> &out,
                               const std::size_t max_num_iter = 50,
                               const double vsat = 100,
                               const Decoder decoder = Decoder::Layered) {
            if (syndrome.size() != n_ra_rows) {
                set_rate(get_n_rows_mother_matrix() - syndrome.size());
            }
            return decode_at_current_rate(llrs, syndrome, out, max_num_iter, vsat, decoder);
        }

        /*!
         * Decode using belief propagation, at the code's current rate (see `decode_infer_rate` to have the
         * rate inferred from `syndrome`'s length instead).
         *
         * @tparam Bit: e.g. std::uint8_t or bool
         * @param llrs: Log likelihood ratios representing the received message
         * @param syndrome: Syndrome of the sent message
         * @param out: Buffer to which the function writes its prediction for the sent message.
         * @param max_num_iter: Maximum number of iterations for the PB algorithm.
         *      Note that the algorithm always terminates automatically when the current prediction matches
         *      the syndrome (early termination), which means that the actual number of iterations cannot be controlled.
         * @param vsat: Cut-off value for messages.
         * @param decoder: which BP variant to run (see `Decoder`); dispatches to `decode_flooding`,
         *      `decode_layered`, or `decode_improved` (gives no access to additional parameters, e.g. damping/rescue settings for `Decoder::Improved`).
         * @return true if and only if the syndrome of buffer `out` matches given `syndrome` (i.e., decoder converged).
         */
        template<std::unsigned_integral Bit>
        bool decode_at_current_rate(const std::vector<double> &llrs,
                                    const std::vector<Bit> &syndrome,
                                    std::vector<Bit> &out,
                                    const std::size_t max_num_iter = 50,
                                    const double vsat = 100,
                                    const Decoder decoder = Decoder::Layered) const {
            switch (decoder) {
                case Decoder::Flooding:
                    return decode_flooding(llrs, syndrome, out, max_num_iter, vsat);
                case Decoder::Layered:
                    return decode_layered(llrs, syndrome, out, max_num_iter, vsat);
                case Decoder::Improved:
                    return decode_improved(llrs, syndrome, out, max_num_iter, vsat);
            }
            throw std::logic_error("decode_at_current_rate: unhandled Decoder");
        }

        /*!
         * Decode using belief propagation with the original FLOODING schedule (all check nodes, then all
         * variable nodes, once per iteration).
         *
         * @tparam Bit: std::uint8_t or bool
         * @param llrs: Log likelihood ratios representing the received message
         * @param syndrome: Syndrome of the sent message
         * @param out: Buffer to which the function writes its prediction for the sent message.
         * @param max_num_iter: Maximum number of iterations for the PB algorithm.
         *      It also terminates when the current prediction matches the syndrome (early termination).
         * @param vsat: Cut-off value for messages.
         * @return whether decoder converged, i.e., whether the syndrome of buffer `out` matches given `syndrome`.
         */
        template<std::unsigned_integral Bit>
        bool decode_flooding(const std::vector<double> &llrs,
                            const std::vector<Bit> &syndrome,
                            std::vector<Bit> &out,
                            const std::size_t max_num_iter = 50,
                            const double vsat = 100) const {
            // check inputs.
            if (llrs.size() != n_cols) {
                throw std::runtime_error("Decoder received invalid input length.");
            }

            if (syndrome.size() != get_n_rows_after_rate_adaption()) {
                throw std::runtime_error(
                        "Decoder (decode_flooding) received invalid syndrome size for current rate. "
                        "Use decode_infer_rate to deduce rate automatically.");
            }

            out.resize(llrs.size());

            std::vector<std::vector<double>> msg_v(n_ra_rows);  // messages from variable nodes to check nodes
            std::vector<std::vector<double>> msg_c(n_cols);  // messages from check nodes to variable nodes

            // initialize msg_v
            for (std::size_t i{}; i < msg_v.size(); ++i) {
                auto &curr_mv = msg_v[i];
                curr_mv.resize(pos_varn[i].size());
                for (std::size_t j{}; j < msg_v[i].size(); ++j) {
                    curr_mv[j] = llrs[pos_varn[i][j]];
                }
            }

            // initialize msg_c
            for (std::size_t i{}; i < msg_c.size(); ++i) {
                msg_c[i].resize(pos_checkn[i].size());
            }

            for (std::size_t it_unused{}; it_unused < max_num_iter; ++it_unused) {
                check_node_update(msg_c, msg_v, syndrome);
                saturate(msg_c, vsat);

                var_node_update(msg_v, msg_c, llrs);
                saturate(msg_v, vsat);

                // hard decision
                hard_decision(out, llrs, msg_c);

                // terminate decoding if codeword matches syndrome
                std::vector<Bit> decision_syndrome(syndrome.size());
                encode_at_current_rate(out, decision_syndrome);
                if (decision_syndrome == syndrome) {
                    return true;
                }

                // check for diverging decoder
                for (const auto &m: msg_v) {
                    for (const auto &v: m) {
                        if (std::isnan(v)) {
                            // TODO maybe use exception?
                            LDPC4QKD_DEBUG_MESSAGE("Decoder Diverged at iteration " << it_unused);
                            return false;
                        }
                    }
                }
            }

            return false;  // Decoding was not successful.
        }

        /*!
         * Decode using belief propagation with LAYERED (serial-C) scheduling.
         * Each check node immediately updates the variable-node posteriors, so information
         * propagates through the graph within a single iteration. Typically converges in
         * roughly half the iterations of flooding and gives equal or better FER.
         */
        template<std::unsigned_integral Bit>
        bool decode_layered(const std::vector<double> &llrs,
                            const std::vector<Bit> &syndrome,
                            std::vector<Bit> &out,
                            const std::size_t max_num_iter = 50,
                            const double vsat = 100) const {
            if (llrs.size() != n_cols) {
                throw std::runtime_error("Decoder received invalid input length.");
            }
            if (syndrome.size() != get_n_rows_after_rate_adaption()) {
                throw std::runtime_error("Decoder received invalid syndrome size for current rate.");
            }
            constexpr double max_tanh = 1. - std::numeric_limits<double>::epsilon();

            out.resize(llrs.size());
            std::vector<double> posterior(llrs);  // current posterior LLR of each variable node

            // check-to-variable message stored per edge, indexed like pos_varn
            std::vector<std::vector<double>> R(n_ra_rows);
            for (std::size_t m{}; m < n_ra_rows; ++m) {
                R[m].assign(pos_varn[m].size(), 0.);
            }

            std::vector<double> tanh_buf;
            std::vector<double> suffix_prod;
            std::vector<Bit> decision_syndrome(syndrome.size());

            for (std::size_t it{}; it < max_num_iter; ++it) {
                for (std::size_t m{}; m < n_ra_rows; ++m) {
                    const auto deg = pos_varn[m].size();
                    tanh_buf.resize(deg);
                    suffix_prod.resize(deg + 1);

                    // variable-to-check messages computed on the fly from current posteriors
                    for (std::size_t k{}; k < deg; ++k) {
                        double t = posterior[pos_varn[m][k]] - R[m][k];
                        t = std::clamp(t, -vsat, vsat);
                        tanh_buf[k] = tanh_half(t);
                    }

                    suffix_prod[deg] = 1.;
                    for (std::size_t k = deg; k-- > 0;) {
                        suffix_prod[k] = suffix_prod[k + 1] * tanh_buf[k];
                    }

                    double prefix = 1 - 2 * static_cast<double>(syndrome[m]);
                    for (std::size_t k{}; k < deg; ++k) {
                        double msg_part = prefix * suffix_prod[k + 1];
                        prefix *= tanh_buf[k];
                        msg_part = std::clamp(msg_part, -max_tanh, max_tanh);
                        const double R_new = std::log1p(msg_part) - std::log1p(-msg_part);

                        // immediately update posterior (this is what makes it "layered")
                        const auto v = pos_varn[m][k];
                        posterior[v] += R_new - R[m][k];
                        R[m][k] = R_new;
                    }
                }

                for (std::size_t j{}; j < n_cols; ++j) {
                    out[j] = posterior[j] < 0 ? 1 : 0;
                }
                encode_at_current_rate(out, decision_syndrome);
                if (decision_syndrome == syndrome) {
                    return true;
                }
            }
            return false;
        }

        /*!
         * Improved decoder: layered SPA with message damping, best-state tracking,
         * and a syndrome-weight bit-flipping rescue stage on failure.
         *
         * \param vsat          cut-off value for messages.
         * \param damping       weight of the new check-to-variable message (1.0 = no damping).
         * \param rescue_weight_cap   only attempt the bit-flip rescue if the best state seen
         *                            has at most this many unsatisfied checks.
         * \param rescue_max_flips    flip budget of the rescue stage.
         */
        template<std::unsigned_integral Bit>
        bool decode_improved(const std::vector<double> &llrs,
                             const std::vector<Bit> &syndrome,
                             std::vector<Bit> &out,
                             const std::size_t max_num_iter = 200,
                             const double vsat = 100,
                             const double damping = 0.8,
                             const std::size_t rescue_weight_cap = 64,
                             const std::size_t rescue_max_flips = 64) const {
            if (llrs.size() != n_cols) {
                throw std::runtime_error("Decoder received invalid input length.");
            }
            if (syndrome.size() != get_n_rows_after_rate_adaption()) {
                throw std::runtime_error("Decoder received invalid syndrome size for current rate.");
            }
            constexpr double max_tanh = 1. - std::numeric_limits<double>::epsilon();

            out.resize(llrs.size());
            std::vector<double> posterior(llrs);
            std::vector<std::vector<double>> R(n_ra_rows);
            for (std::size_t m{}; m < n_ra_rows; ++m) {
                R[m].assign(pos_varn[m].size(), 0.);
            }

            std::vector<double> tanh_buf, suffix_prod;
            std::vector<Bit> decision_syndrome(syndrome.size());
            std::vector<Bit> best_out;
            std::size_t best_weight = std::numeric_limits<std::size_t>::max();

            for (std::size_t it{}; it < max_num_iter; ++it) {
                for (std::size_t m{}; m < n_ra_rows; ++m) {
                    const auto deg = pos_varn[m].size();
                    tanh_buf.resize(deg);
                    suffix_prod.resize(deg + 1);
                    for (std::size_t k{}; k < deg; ++k) {
                        double t = posterior[pos_varn[m][k]] - R[m][k];
                        t = std::clamp(t, -vsat, vsat);
                        tanh_buf[k] = tanh_half(t);
                    }
                    suffix_prod[deg] = 1.;
                    for (std::size_t k = deg; k-- > 0;) {
                        suffix_prod[k] = suffix_prod[k + 1] * tanh_buf[k];
                    }
                    double prefix = 1 - 2 * static_cast<double>(syndrome[m]);
                    for (std::size_t k{}; k < deg; ++k) {
                        double msg_part = prefix * suffix_prod[k + 1];
                        prefix *= tanh_buf[k];
                        msg_part = std::clamp(msg_part, -max_tanh, max_tanh);
                        const double R_bp = std::log1p(msg_part) - std::log1p(-msg_part);
                        // damping: convex combination of old and new message
                        const double R_new = damping * R_bp + (1. - damping) * R[m][k];
                        const auto v = pos_varn[m][k];
                        posterior[v] += R_new - R[m][k];
                        R[m][k] = R_new;
                    }
                }

                for (std::size_t j{}; j < n_cols; ++j) {
                    out[j] = posterior[j] < 0 ? 1 : 0;
                }
                encode_at_current_rate(out, decision_syndrome);

                // best-state tracking: count unsatisfied checks
                std::size_t weight = 0;
                for (std::size_t m{}; m < decision_syndrome.size(); ++m) {
                    weight += (decision_syndrome[m] != syndrome[m]);
                }
                if (weight == 0) {
                    return true;
                }
                if (weight < best_weight) {
                    best_weight = weight;
                    best_out = out;
                }
            }

            LDPC4QKD_DEBUG_MESSAGE("[decode_improved] best_weight at failure: " << best_weight);
            // ---- bit-flip rescue stage on the best state seen ----
            if (best_weight > rescue_weight_cap) {
                out = best_out.empty() ? out : best_out;
                return false;
            }
            out = best_out;
            encode_at_current_rate(out, decision_syndrome);

            // column -> adjacent (rate-adapted) check nodes, built once per rescue
            std::vector<std::vector<idx_t>> col_to_checks(n_cols);
            for (std::size_t m{}; m < n_ra_rows; ++m) {
                for (const auto v : pos_varn[m]) {
                    col_to_checks[v].push_back(static_cast<idx_t>(m));
                }
            }

            std::vector<std::uint8_t> unsat(n_ra_rows);
            std::size_t weight = 0;
            for (std::size_t m{}; m < n_ra_rows; ++m) {
                unsat[m] = (decision_syndrome[m] != syndrome[m]);
                weight += unsat[m];
            }

            std::vector<std::uint8_t> flipped(n_cols, 0);  // tabu marker for escape moves
            for (std::size_t flip{}; flip < rescue_max_flips && weight > 0; ++flip) {
                // candidates: variables adjacent to at least one unsatisfied check.
                // Greedy: take the best positive-gain flip. If none exists (absorbing-set
                // pattern: every wrong bit sees a majority of satisfied checks), take an
                // escape move: the not-yet-flipped candidate with the lowest reliability.
                long best_gain = std::numeric_limits<long>::min();
                double best_rel = std::numeric_limits<double>::infinity();
                std::size_t best_v = n_cols;
                long esc_gain = std::numeric_limits<long>::min();
                double esc_rel = std::numeric_limits<double>::infinity();
                std::size_t esc_v = n_cols;
                for (std::size_t m{}; m < n_ra_rows; ++m) {
                    if (!unsat[m]) continue;
                    for (const auto v : pos_varn[m]) {
                        long u = 0;
                        for (const auto c : col_to_checks[v]) u += unsat[c];
                        const long gain = 2 * u - static_cast<long>(col_to_checks[v].size());
                        const double rel = std::abs(posterior[v]);
                        if (gain > best_gain || (gain == best_gain && rel < best_rel)) {
                            best_gain = gain;
                            best_rel = rel;
                            best_v = v;
                        }
                        if (!flipped[v] &&
                            (rel < esc_rel || (rel == esc_rel && gain > esc_gain))) {
                            esc_gain = gain;
                            esc_rel = rel;
                            esc_v = v;
                        }
                    }
                }
                std::size_t v_flip;
                if (best_v != n_cols && best_gain > 0) {
                    v_flip = best_v;            // strict descent
                } else if (esc_v != n_cols) {
                    v_flip = esc_v;             // sideways/uphill escape, tabu-guarded
                } else {
                    break;                      // nothing left to try
                }
                flipped[v_flip] = 1;
                out[v_flip] = out[v_flip] ? 0 : 1;
                for (const auto c : col_to_checks[v_flip]) {
                    if (unsat[c]) { weight--; } else { weight++; }
                    unsat[c] = !unsat[c];
                }
            }
            return weight == 0;
        }

        //! manually trigger rate adaption. In normal circumstances, the user does not need this function
        //! \param n_line_combs number of line combinations to use (starting from the mother code)
        void set_rate(std::size_t n_line_combs) {
            recompute_pos_vn_cn(n_line_combs);
        }

        template<typename BitL=bool, typename BitR=bool>
        constexpr void encode_at_current_rate(
                const std::vector<BitL> &in, std::vector<BitR> &out) const {
            if (in.size() != n_cols) {
                LDPC4QKD_DEBUG_MESSAGE("Encoder received invalid input length.");  // TODO maybe use exception?
                return;
            }

            out.assign(pos_varn.size(), 0);

            for (std::size_t i{}; i < pos_varn.size(); ++i) {
                for (auto &var_node: pos_varn[i]) {
                    out[i] = xor_as_bools(out[i], in[var_node]);
                }
            }
        }

        bool operator==(const RateAdaptiveCode &rhs) const {
            return n_mother_rows == rhs.n_mother_rows &&
                   n_cols == rhs.n_cols &&
                   mother_pos_varn == rhs.mother_pos_varn &&
                   rows_to_combine == rhs.rows_to_combine &&
                   pos_checkn == rhs.pos_checkn &&
                   pos_varn == rhs.pos_varn &&
                   n_ra_rows == rhs.n_ra_rows;
        }

        // ----------------------------------------------------------------------------------------- getters and setters
        [[nodiscard]]
        const std::vector<std::vector<idx_t>> &getPosCheckn() const {
            return pos_checkn;
        }

        [[nodiscard]]
        const std::vector<std::vector<idx_t>> &getPosVarn() const {
            return pos_varn;
        }

        /// ignores rate adaption! Only gives number of rows in the mother matrix.
        [[nodiscard]] auto get_n_rows_mother_matrix() const {
            return n_mother_rows;
        }

        /// Includes rate adaption. Access to internal state!
        [[nodiscard]] auto get_n_rows_after_rate_adaption() const {
            return n_ra_rows;
        }

        [[nodiscard]] auto getNCols() const {
            return n_cols;
        }

        //! Maximum number of line combinations available for rate adaption, i.e. the largest `n_line_combs` that
        //! `set_rate()` accepts, counted as steps of one combined pair each starting from the mother matrix (0
        //! line combs = the unmodified mother matrix). Either from the auto-generated (LCG-based) scheme (used
        //! when `rows_to_combine` is empty) or from the explicit `rows_to_combine`.
        [[nodiscard]] std::size_t get_max_ra_steps() const {
            return rows_to_combine.empty() ? (n_mother_rows / 2) : (rows_to_combine.size() / 2);
        }

    private:   // -------------------------------------------------------------------------------------- private members
        template<typename BitL, typename BitR>
        constexpr static bool xor_as_bools(BitL lhs, BitR rhs) {
            return (static_cast<bool>(lhs) != static_cast<bool>(rhs));
        }

        template<typename Idx>
        static Idx compute_n_cols(std::vector<std::vector<Idx>> mother_pos_varn) {
            if (mother_pos_varn.empty()) {
                return 0;
            } else {
                Idx result{};
                for (const auto &v: mother_pos_varn) {
                    auto current_max = *std::max_element(v.cbegin(), v.cend());
                    result = std::max(result, current_max);
                }
                return result + 1; // add one because indices in `mother_pos_varn` are zero-based.
            }
        }

        /*!
         * compute `mother_pos_varn` from `colptr` and `rowIdx` for a given LDPC matrix stored in compressed sparse column
         * format. "Values" array is omitted because all values are assumed to be 1 (binary LDPC matrix).
         *
         * @tparam idx_t unsigned integer type fitting number of columns N (thus also number of rows M)
         * @tparam colptr_t unsigned integer type that fits ("number of non-zero matrix entries" + 1)
         * @param colptr column pointer array for specifying mother parity check matrix.
         * @param rowIdx row index array for specifying mother parity check matrix.
         * @return Input variable nodes to each check node (of the Tanner graph)
         */
        template<typename colptr_t>
        static std::vector<std::vector<idx_t>> compute_mother_pos_varn(
                const std::vector<colptr_t> &colptr,
                const std::vector<idx_t> &rowIdx) {
            // number of columns in full matrix represented by given compressed sparse column (CSC) storage
            const auto nCols = colptr.size() - 1;
            // number of rows in full matrix represented by given compressed sparse column (CSC) storage
            const auto nMotherRows = *std::max_element(rowIdx.begin(), rowIdx.end()) + 1u;

            std::vector<std::vector<idx_t>> pos_varn_tmp{nMotherRows, std::vector<idx_t>{}};
            for (idx_t col = 0; col < nCols; col++) {
                for (auto j = colptr[col]; j < colptr[col + 1u]; j++) {
                    pos_varn_tmp[rowIdx[j]].push_back(col);
                }
            }
            return pos_varn_tmp;
        }

        template<typename Bit>
        void check_node_update(std::vector<std::vector<double>> &msg_c,
                               const std::vector<std::vector<double>> &msg_v,
                               const std::vector<Bit> &syndrome) const {
            // Largest value strictly below 1, so that log1p(x) - log1p(-x) stays finite.
            constexpr double max_tanh = 1. - std::numeric_limits<double>::epsilon();
            std::vector<idx_t> mc_position(n_cols);
            std::vector<double> tanh_buf;    // tanh(msg/2) of each incoming message
            std::vector<double> suffix_prod; // suffix products of tanh_buf

            for (std::size_t m{}; m < n_ra_rows; ++m) {
                // Note: pos_varn[m].size() = check_node_degrees[m]
                const auto curr_check_node_degree = pos_varn[m].size();
                tanh_buf.resize(curr_check_node_degree);
                suffix_prod.resize(curr_check_node_degree + 1);

                for (std::size_t k{}; k < curr_check_node_degree; ++k) {
                    tanh_buf[k] = tanh_half(msg_v[m][k]);
                }

                // suffix_prod[k] = product of tanh_buf[k..deg-1]
                suffix_prod[curr_check_node_degree] = 1.;
                for (std::size_t k = curr_check_node_degree; k-- > 0;) {
                    suffix_prod[k] = suffix_prod[k + 1] * tanh_buf[k];
                }

                // prefix accumulates syndrome sign times product of tanh_buf[0..k-1]
                double prefix = 1 - 2 * static_cast<double>(syndrome[m]);
                for (std::size_t k{}; k < curr_check_node_degree; ++k) {
                    // extrinsic product: all incoming tanh's except edge k. No division needed,
                    // so exact zeros and saturated (+-1) messages are handled correctly.
                    double msg_part = prefix * suffix_prod[k + 1];
                    prefix *= tanh_buf[k];

                    // clamp into (-1, 1) so the result is always finite (cf. AFF3CT SPA decoder)
                    msg_part = std::clamp(msg_part, -max_tanh, max_tanh);

                    // log1p is more accurate than log((1+x)/(1-x)) for |msg_part| << 1
                    const double msg_final = std::log1p(msg_part) - std::log1p(-msg_part);

                    // place the message at the correct position in the output array
                    const idx_t curr_pos_varn = pos_varn[m][k];
                    msg_c[curr_pos_varn][mc_position[curr_pos_varn]] = msg_final;
                    mc_position[curr_pos_varn]++;
                }
            }
        }

        void var_node_update(std::vector<std::vector<double>> &msg_v,
                             const std::vector<std::vector<double>> &msg_c,
                             const std::vector<double> &llrs) const {
            std::vector<idx_t> mv_position(n_cols);

            for (std::size_t m{}; m < llrs.size(); ++m) {
                const double mv_sum = std::accumulate(msg_c[m].begin(), msg_c[m].end(), llrs[m]);

                // Note: pos_checkn[m].size() = var_node_degs[m]
                for (std::size_t k{}; k < pos_checkn[m].size(); ++k) {
                    const double msg = mv_sum - msg_c[m][k];

                    // place the message at the correct position in the output array
                    const idx_t curr_pos_cn = pos_checkn[m][k];
                    msg_v[curr_pos_cn][mv_position[curr_pos_cn]] = msg;
                    mv_position[curr_pos_cn]++;
                }
            }
        }

        template<typename Bit=bool>
        void hard_decision(
                std::vector<Bit> &out,
                const std::vector<double> &llrs,
                const std::vector<std::vector<double>> &msg_c) const {
            std::fill(out.begin(), out.end(), 0);
            for (std::size_t j{}; j < llrs.size(); ++j) {
                const double curr_sum = std::accumulate(msg_c[j].begin(), msg_c[j].end(), llrs[j]);
                if (curr_sum < 0) {
                    out[j] = 1;
                }
            }
        }

        template<typename T>
        static void saturate(std::vector<std::vector<T>> &mv, const T vsat) {
            for (auto &v: mv) {
                for (auto &a: v) {
                    if (a > vsat) { a = vsat; }
                    else if (a < -vsat) { a = -vsat; }
                }
            }
        }

        /*!
         * Recompute inner representation of rate adapted LDPC code (`pos_varn` and `pos_cn`),
         * starting from the mother code represented by `mother_pos_varn`.
         *
         * The combined row's variable-node list is built as the union of the two input rows'
         * lists. This corresponds to elementwise OR of the two rows.
         * However, this is correct for BOTH row-selection methods:
         *   - auto-generated (LCG-based, `rows_to_combine.empty()`): pairs are not guaranteed to
         *     be variable-disjoint (collisions do occur in practice, see `encode_with_ra`'s own
         *     comment), so OR is the only choice that never silently drops a shared variable's
         *     connection to the Tanner graph (XOR/symmetric-difference would, whenever both rows
         *     already touch that variable).
         *   - explicit (`rows_to_combine` given): these pairs are meant to be pre-verified
         *     variable-disjoint, so union and symmetric-difference always agree.
         * See `encode_with_ra` for the corresponding syndrome-*value* combination (which, unlike
         * this structural step, does need to differ by scheme: XOR for explicit, OR for
         * auto-generated).
         *
         * @param n_line_combs number of line combinations to perform for rate adaption.
         */
        void recompute_pos_vn_cn(std::size_t n_line_combs) {
            if (get_max_ra_steps() < n_line_combs) {
                throw std::runtime_error("Requested rate not supported. Not enough line combinations specified.");
            }

            {   // recompute pos_varn ---------------------------------------------------------------
                // This uses different size vectors for nodes with different degrees.
                // Alternatively, one could set the sizes to be the same (set them to the largest check node degree)

                n_ra_rows = n_mother_rows - n_line_combs;
                pos_varn.assign(n_ra_rows, std::vector<idx_t>{});

                if (n_line_combs == 0) {
                    pos_varn = mother_pos_varn;
                } else {
                    // Make temporary copy of `mother_pos_varn`
                    std::vector<std::vector<idx_t>> pos_varn_nora{mother_pos_varn};

                    // If `rows_to_combine` is empty, pairs are auto-generated (LCG-based) and placed at the FRONT
                    // of the output; otherwise the explicitly given pairs are used, placed at the BACK.
                    const bool auto_generated = rows_to_combine.empty();
                    RateAdaptLCG::LCG local_lcg = auto_rate_adaption;

                    // for auto-generated rate adaption, combined lines go at the front: [0, n_line_combs).
                    // for explicit rate adaption, combined lines go at the back: [start_of_ra_part, n_ra_rows).
                    const auto start_of_ra_part = auto_generated ? 0 : (n_mother_rows - 2 * n_line_combs);

                    for (std::size_t i{}; i < n_line_combs; ++i) {
                        idx_t idx1, idx2;
                        if (auto_generated) {
                            idx1 = static_cast<idx_t>(local_lcg.next());
                            idx2 = static_cast<idx_t>(local_lcg.next());
                        } else {
                            idx1 = rows_to_combine[2 * i];
                            idx2 = rows_to_combine[2 * i + 1];
                        }

                        auto &curr_varn_vec = pos_varn[start_of_ra_part + i];
                        curr_varn_vec.insert(curr_varn_vec.end(),
                                             pos_varn_nora[idx1].begin(),
                                             pos_varn_nora[idx1].end());
                        curr_varn_vec.insert(curr_varn_vec.end(),
                                             pos_varn_nora[idx2].begin(),
                                             pos_varn_nora[idx2].end());

                        pos_varn_nora[idx1].clear();
                        pos_varn_nora[idx2].clear();

                        // Union of the two rows' variable-node lists, i.e., elementwise OR of rows
                        // Correct for both schemes, see docstring.
                        std::sort(curr_varn_vec.begin(), curr_varn_vec.end());
                        curr_varn_vec.erase(std::unique(curr_varn_vec.begin(), curr_varn_vec.end()),
                                            curr_varn_vec.end());
                    }

                    std::size_t j{};

                    // put the remaining lines that were not rate adapted into the rest of the new LDPC code
                    // (after the combined lines for auto-generated rate adaption; before them otherwise).
                    const std::size_t leftover_begin = auto_generated ? n_line_combs : 0;
                    const std::size_t leftover_end = auto_generated ? n_ra_rows : start_of_ra_part;
                    for (std::size_t i = leftover_begin; i < leftover_end; ++i) {
                        while (pos_varn_nora[j].empty()) {
                            j++;
                        }
                        pos_varn[i] = std::move(pos_varn_nora[j]);
                        j++;
                    }

                }
            }  // end recompute pos_varn

            {   // recompute pos_checkn -------------------------------------------------------------------------------
                // Now compute pos_checkn from the previously computed pos_varn.
                // These arrays contain the same information.
                pos_checkn.assign(n_cols, std::vector<idx_t>{});

                for (idx_t i{}; i < pos_varn.size(); ++i) {
                    for (auto &vn: pos_varn[i]) {
                        pos_checkn[vn].push_back(i);
                    }
                }
            }  // end recompute pos_checkn
        }

        // ---------------------------------------------------------------------------------------------- private fields
        // TODO consider making these of type `idx_t`.
        const std::size_t n_mother_rows;  // const because it's not possible to change the mother matrix
        const std::size_t n_cols;  // const because it's not possible to change the mother matrix

        /// Input variable nodes to each check node of the mother matrix.
        /// Rate adaption always starts from here.
        /// Can be obtained from
        /// (1) arrays `colptr` and `row_idx`, which represent the binary LDPC matrix in CSC format
        /// (2) an "encoder" implementing `ComputablePosVar`.
        const std::vector<std::vector<idx_t>> mother_pos_varn;

        /// stores specification of rate adaption.
        /// Each rate adaption is re-computed using `mother_pos_checkn` and `rows_to_combine`.
        /// If empty, rate adaption pairs are instead generated on demand from `auto_rate_adaption`.
        const std::vector<idx_t> rows_to_combine;

        /// Initial state of the LCG used to auto-generate rate adaption pairs whenever `rows_to_combine` is empty.
        /// A pure deterministic function of `n_mother_rows` (fixed seed 0), computed unconditionally, but only
        /// actually consumed (by replaying it from this state) when `rows_to_combine` is empty.
        const RateAdaptLCG::LCG auto_rate_adaption;

        /// `pos_checkn` and `pos_varn` store the current rate adapted code, which is actually used for decoding.
        std::vector<std::vector<idx_t>> pos_checkn;  /// Input check nodes to each variable node
        std::vector<std::vector<idx_t>> pos_varn;  /// Input variable nodes to each check node

        /// current number of matrix rows (given current rate adaption).
        std::size_t n_ra_rows{};
    };

}

#endif //LDPC4QKD_LDPC_MATRIX_HPP
