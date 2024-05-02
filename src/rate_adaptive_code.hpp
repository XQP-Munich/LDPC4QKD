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
#include <exception>
#include <stdexcept>

#ifdef LDPC4QKD_DEBUG_MESSAGES_ENABLED

#include <iostream>

#define LDPC4QKD_DEBUG_MESSAGE(msg) do {std::cerr << msg << std::endl;} while (false)

#else

#define LDPC4QKD_DEBUG_MESSAGE(msg)

#endif /* ifdef LDPC4QKD_DEBUG_MESSAGES_ENABLED */


namespace LDPC4QKD {

    /*!
     * compute log-likelihood-ratios for a given keye and channel parameter
     * @tparam Bit e.g. std::uint8_t or bool
     * @param bitstring string of bits
     * @param bsc_channel_parameter bit-flip probability of binary symmetric channel (BSC)
     * @return log-likelyhoods corresponding to input bitstring
     */
    template<typename Bit>
    std::vector<double> llrs_bsc(const std::vector<Bit> &bitstring, const double bsc_channel_parameter) {
        double vlog = log((1 - bsc_channel_parameter) / bsc_channel_parameter);
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
     * (TODO use concept `std::unsigned_integral` when using C++20)
     *
     * @tparam colptr_t unsigned integert type that fits ("number of non-zero matrix entries" + 1)
     * @tparam idx_t unsigned integer type fitting number of columns N (thus also number of rows M)
     */
    template<typename colptr_t=std::uint32_t,
            typename idx_t=std::uint16_t>
    class RateAdaptiveCode {
    public:
        // ------------------------------------------------------------------------------------------------ type aliases
        using MatrixIndex = idx_t;
        using ColumnPointer = colptr_t;

        // ------------------------------------------------------------------------------------------------ constructors
        /*!
         * Constructor for using the code without rate adaption.
         * The parity check matrix matrix is stored using Compressed Sparse Column (CSC) format.
         *
         * @param colptr column pointer array for specifying mother parity check matrix.
         * @param rowIdx row index array for specifying mother parity check matrix.
         */
        RateAdaptiveCode(const std::vector<colptr_t> &colptr, const std::vector<idx_t> &rowIdx)
                : n_mother_rows(*std::max_element(rowIdx.begin(), rowIdx.end()) + 1u),
                  n_cols(colptr.size() - 1),
                  mother_pos_varn(compute_mother_pos_varn(colptr, rowIdx)), // computed here and henceforth `const`!
                  rows_to_combine({}) {
            constexpr idx_t n_line_combs = 0;
            recompute_pos_vn_cn(n_line_combs);
        }

        /*!
         * Constructor for using the code with rate adaption.
         * The mother parity check matrix is stored using Compressed Sparse Column (CSC) format.
         * The rate adaption is stored as an array of matrix row indices, which are combined for rate adaption.
         *
         * note: if you for some reason find yourself creating a lot of such objects and you know that there are no
         * variable node eliminations, disable the elimination check to speed up this constructor.
         *
         * note: invalid `rows_to_combine_rate_adapt`, for example non-zero based, may lead to a segmentation fault.
         * TODO add extra checks for validity of `rows_to_combine_rate_adapt`
         *
         * note: there used to be a parameter `do_elimination_check` to check for repeated node indices after rate adaption.
         *      Such indices are now removed during `recompute_pos_vn_cn`. Consequentially, node eliminations are allowed.
         *              TODO reconsider this and remove commented-out function `has_var_node_eliminations` below
         *
         * @param colptr column pointer array for specifying mother parity check matrix.
         * @param rowIdx row index array for specifying mother parity check matrix.
         * @param rows_to_combine_rate_adapt array of mother-matrix line indices to be combined for rate adaption
         * @param initial_row_combs number of line indices to combine initially
         */
        RateAdaptiveCode(std::vector<colptr_t> colptr,
                         std::vector<idx_t> rowIdx,
                         std::vector<idx_t> rows_to_combine_rate_adapt,
                         idx_t initial_row_combs = 0)
                : n_mother_rows(*std::max_element(rowIdx.begin(), rowIdx.end()) + 1u),
                  n_cols(colptr.size() - 1),
                  mother_pos_varn(compute_mother_pos_varn(colptr, rowIdx)), // computed here and henceforth `const`!
                  rows_to_combine(std::move(rows_to_combine_rate_adapt)) {
            if (rows_to_combine.size() % 2 != 0) {
                throw std::domain_error("The number of rows to combine for rate adaption "
                                        "(size of argument array) is an odd number (expected even).");
            }

            if (initial_row_combs > rows_to_combine.size() / 2) {
                throw std::domain_error("The number of desired initial row combinations for rate adaption "
                                        "is larger than the given array of lines to combine.");
            }

            // compute current `pos_varn` and `pos_checkn` from `mother_pos_varn`
            recompute_pos_vn_cn(initial_row_combs);

//            if (do_elimination_check && has_var_node_eliminations()) {
//                throw std::domain_error("Given rate adaption implies variable node eliminations. "
//                                        "Rate adaption with eliminations degrades performance. Do not use!");
//            }
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
                  n_ra_rows(n_mother_rows - initial_row_combs) {
            if (rows_to_combine.size() % 2 != 0) {
                throw std::domain_error("The number of rows to combine for rate adaption "
                                        "(size of argument array) is an odd number (expected even).");
            }

            if (initial_row_combs > rows_to_combine.size() / 2) {
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
         * @tparam Bit e.g. std::uint8_t or bool. TODO use concept `std::unsigned_integral` when using C++20
         * @param in input array
         * @param out Vector to store syndrome. Will be resized to `output_syndrome_length`
         * @param output_syndrome_length Desired length of syndrome (exception is thrown if not satisfiable)
         */
        template<typename Bit>
        void encode_with_ra(
                const std::vector<Bit> &in, std::vector<Bit> &out, std::size_t output_syndrome_length) const {
            if (in.size() != n_cols) {
                throw std::domain_error("Encoder (encode_with_ra) received invalid input length.");
            }
            if (output_syndrome_length > n_mother_rows) {
                throw std::domain_error("Requested syndrome is larger than the number of rows of the mother matrix.");
            }
            if (output_syndrome_length < n_mother_rows - (rows_to_combine.size() / 2)) {
                throw std::domain_error("Requested syndrome is smaller than supported by the specified rate adaption.");
            }

            // HAVE to use !!!SIGNED!!! `int8_t`!!! Value `-1` is used below to mark bits included into final output!
            std::vector<std::int8_t> non_ra_encoding;
            encode_no_ra(in, non_ra_encoding);

            // now use the non-rate adapted syndrome to compute the rate adapted syndrome
            const std::size_t n_line_combinations = n_mother_rows - output_syndrome_length;
            out.assign(output_syndrome_length, 0);

            std::size_t start_of_ra_part = output_syndrome_length - n_line_combinations;

            // put results of combined lines at the back of output.
            for (std::size_t i{}; i < n_line_combinations; ++i) {
                out[start_of_ra_part + i] = xor_as_bools(non_ra_encoding[rows_to_combine[2 * i]],
                                                         non_ra_encoding[rows_to_combine[2 * i + 1]]);
                non_ra_encoding[rows_to_combine[2 * i]] = -1;  // -1 marks that the value has been used.
                non_ra_encoding[rows_to_combine[2 * i + 1]] = -1;
            }

            std::size_t j{};
            // put the remaining bits that were not rate adapted at the front of output.
            for (std::size_t i{}; i < start_of_ra_part; ++i) {
                while (non_ra_encoding[j] == -1) {
                    j++;
                }
                out[i] = non_ra_encoding[j];
                j++;
            }
        }

        /// decoder infers rate from the length of the syndrome and changes the internal decoder state to match this rate.
        /// Note: since this function modifies the code (by performing rate adaption), it is NOT CONST.
        /// this change may be somewhat computationally expensive TODO benchmark this
        /// `Bit` should be e.g. std::uint8_t or bool. TODO use concept `std::unsigned_integral` when using C++20
        template<typename Bit>
        bool decode_infer_rate(const std::vector<double> &llrs,
                               const std::vector<Bit> &syndrome,
                               std::vector<Bit> &out,
                               const std::size_t max_num_iter = 50,
                               const double vsat = 100) {
            if (syndrome.size() != n_ra_rows) {
                set_rate(get_n_rows_mother_matrix() - syndrome.size());
            }
            return decode_at_current_rate(llrs, syndrome, out, max_num_iter, vsat);
        }


        /*!
         * Decode using belief propagation
         *
         * @tparam Bit: e.g. std::uint8_t or bool TODO use concept `std::unsigned_integral` when using C++20
         * @param llrs: Log likelihood ratios representing the received message
         * @param syndrome: Syndrome of the sent message
         * @param out: Buffer to which the function writes its prediction for the sent message.
         * @param max_num_iter: Maximum number of iterations for the PB algorithm.
         *      Note that the algorithm always terminates automatically when the current prediction matches
         *      the syndrome (early termination), which means that the actual number of iterations cannot be controlled.
         * @param vsat: Cut-off value for messages.
         * @return true if and only if the syndrome of buffer `out` matches given `syndrome` (i.e., decoder converged).
         */
        template<typename Bit>
        bool decode_at_current_rate(const std::vector<double> &llrs,
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
                        "Decoder (decode_at_current_rate) received invalid syndrome size for current rate. "
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

        [[nodiscard]] auto get_max_ra_steps() const {
            return rows_to_combine.size() / 2;
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

        /// compute `mother_pos_varn` from `colptr` and `rowIdx`
        static std::vector<std::vector<idx_t>> compute_mother_pos_varn(
                const std::vector<colptr_t> &colptr,
                const std::vector<idx_t> &rowIdx) {
            // number of columns in full matrix represented by given compressed sparse column (CSC) storage
            const auto n_cols = colptr.size() - 1;
            // number of rows in full matrix represented by given compressed sparse column (CSC) storage
            const auto n_mother_rows = *std::max_element(rowIdx.begin(), rowIdx.end()) + 1u;

            std::vector<std::vector<idx_t>> pos_varn_tmp{n_mother_rows, std::vector<idx_t>{}};
            for (idx_t col = 0; col < n_cols; col++) {
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
            double msg_part{};
            std::vector<idx_t> mc_position(n_cols);

            for (std::size_t m{}; m < n_ra_rows; ++m) {
                // product of incoming messages
                double mc_prod = 1 - 2 * static_cast<double>(syndrome[m]);
                // Note: pos_varn[m].size() = check_node_degrees[m]
                const auto curr_check_node_degree = pos_varn[m].size();
                for (std::size_t k{}; k < curr_check_node_degree; ++k) {
                    mc_prod *= ::tanh(0.5 * msg_v[m][k]);
                }

                for (std::size_t k{}; k < curr_check_node_degree; ++k) {
                    // computing message from
                    if (msg_v[m][k] == 0.) {  // TODO test this bit more carefully.
                        LDPC4QKD_DEBUG_MESSAGE("Decoder went into the 'untested bit'!!");
                        msg_part = 1;
                        for (std::size_t non_k{}; non_k < curr_check_node_degree; ++non_k) {
                            if (non_k != k) {
                                msg_part *= ::tanh(0.5 * msg_v[m][k]);
                            }
                        }
                    } else {
                        msg_part = mc_prod / ::tanh(0.5 * msg_v[m][k]);
                    }

                    auto msg_final = ::log((1 + msg_part) / (1 - msg_part));

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
        static void saturate(std::vector<std::vector<T>> &mv, const double vsat) {
            for (auto &v : mv) {
                for (auto &a : v) {
                    if (a > vsat) { a = vsat; }
                    else if (a < -vsat) { a = -vsat; }
                }
            }
        }

        /*!
         * Recompute inner representation of rate adapted LDPC code (`pos_varn` and `pos_cn`),
         * starting from the mother code represented by `mother_pos_varn`.
         * Note: this function "deals incorrectly" with variable node elimination during rate adaption.
         * variable node elimination should not happen in the first place
         *
         * @param n_line_combs number of line combinations to perform for rate adaption.
         */
        void recompute_pos_vn_cn(std::size_t n_line_combs) {
            if (rows_to_combine.size() < 2 * n_line_combs) {
                throw std::runtime_error("Requested rate not supported. Not enough line combinations specified.");
            }

            {   // recompute pos_varn ---------------------------------------------------------------
                // TODO check if this assumes full rank of H (should have that anyway)
                // This uses different size vectors for nodes with different degrees.
                // Alternatively, one could set the sizes to be the same (set them to the largest check node degree)

                n_ra_rows = n_mother_rows - n_line_combs;
                pos_varn.assign(n_ra_rows, std::vector<idx_t>{});

                if (n_line_combs == 0) {
                    pos_varn = mother_pos_varn;
                } else {
                    // Make temporary copy of `mother_pos_varn`
                    std::vector<std::vector<idx_t>> pos_varn_nora{mother_pos_varn};

                    // put results of combined lines at the back of the new LDPC code
                    const auto start_of_ra_part = n_mother_rows - 2 * n_line_combs;

                    for (std::size_t i{}; i < n_line_combs; ++i) {
                        auto &curr_varn_vec = pos_varn[start_of_ra_part + i];
                        curr_varn_vec.insert(curr_varn_vec.end(),
                                             pos_varn_nora[rows_to_combine[2 * i]].begin(),
                                             pos_varn_nora[rows_to_combine[2 * i]].end());
                        curr_varn_vec.insert(curr_varn_vec.end(),
                                             pos_varn_nora[rows_to_combine[2 * i + 1]].begin(),
                                             pos_varn_nora[rows_to_combine[2 * i + 1]].end());

                        pos_varn_nora[rows_to_combine[2 * i]].clear();
                        pos_varn_nora[rows_to_combine[2 * i + 1]].clear();

                        // TODO speed up this part by producing the rate adaption as unique positions and already sorted
                        std::sort(curr_varn_vec.begin(), curr_varn_vec.end());
                        curr_varn_vec.erase(std::unique(curr_varn_vec.begin(), curr_varn_vec.end()),
                                            curr_varn_vec.end());
                    }

                    std::size_t j{};

                    // put the remaining lines that were not rate adapted at the front of the new LDPC code.
                    for (std::size_t i{}; i < start_of_ra_part; ++i) {
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
                    for (auto &vn : pos_varn[i]) {
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
        const std::vector<idx_t> rows_to_combine;  // TODO if empty, use random rate adaption

        /// `pos_checkn` and `pos_varn` store the current rate adapted code, which is actually used for decoding.
        std::vector<std::vector<idx_t>> pos_checkn;  /// Input check nodes to each variable node
        std::vector<std::vector<idx_t>> pos_varn;  /// Input variable nodes to each check node

        /// current number of matrix rows (given current rate adaption).
        std::size_t n_ra_rows{};
    };

}

#endif //LDPC4QKD_LDPC_MATRIX_HPP
