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
#include <set>

#ifdef DEBUG_MESSAGES_ENABLED

#include <iostream>

#define DEBUG_MESSAGE(msg) do {std::cerr << msg << std::endl;} while (false)

#else

#define DEBUG_MESSAGE(msg)

#endif /* ifdef DEBUG_MESSAGES_ENABLED */


namespace LDPC4QKD {

    /*!
     * Belief propagation (BP) decoder for binary low density parity check (LDPC) codes.
     * Supports rate adaption (reducing the number of LDPC matrix rows).
     * Intended for distributed source coding (a.k.a. Slepian-Wolf coding).
     * LDPC code is stored in sparse column storage (CSC) format.
     *
     * @tparam Bit type of matrix entires (only values zero and one are used, default to bool)
     * @tparam colptr_t unsigned integer type that fits ("number of non-zero matrix entries" + 1)
     * @tparam idx_t unsigned integer type fitting number of columns N (thus also number of rows M)
     */
    template<typename Bit,
            typename colptr_t=std::uint32_t,
            typename idx_t=std::uint16_t>
    class RateAdaptiveCode {
    public:
        // ------------------------------------------------------------------------------------------------ type aliases
        using NodeDegree = std::uint16_t;
        using MatrixIndex = idx_t;
        using ColumnPointer = colptr_t;
        using MatrixEntry = Bit;

        static_assert(!std::numeric_limits<NodeDegree>::is_signed, "Expecting unsigned integer types");
        static_assert(!std::numeric_limits<MatrixIndex>::is_signed, "Expecting unsigned integer types");
        static_assert(!std::numeric_limits<ColumnPointer>::is_signed, "Expecting unsigned integer types");
        static_assert(std::numeric_limits<MatrixEntry>::is_integer,
                "Expecting integer type (e.g. bool or uint8) for the bit values. Will use only 'true' and 'false'.");

        // ------------------------------------------------------------------------------------------------ constructors
        /*!
         * Constructor for using the code without rate adaption.
         * The parity check matrix matrix is stored using Compressed Sparse Column (CSC) format.
         *
         * @param colptr column pointer array for specifying mother parity check matrix.
         * @param rowIdx row index array for specifying mother parity check matrix.
         */
        RateAdaptiveCode(const std::vector<colptr_t> &colptr, const std::vector<idx_t> &rowIdx)
                : colptr(colptr), row_idx(rowIdx), rows_to_combine({}),
                  n_mother_rows(*std::max_element(rowIdx.begin(), rowIdx.end()) + 1u),
                  n_cols(colptr.size() - 1),
                  n_ra_rows(n_mother_rows) {
            constexpr std::size_t n_line_combs = 0;

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
         *      Such indices are now removed during `recompute_pos_vn_cn`. Consequentially, node elliminations are allowed.
         *              TODO reconsider this and remove commented-out function `has_var_node_eliminations` below
         *
         * @param colptr colptr column pointer array for specifying mother parity check matrix.
         * @param rowIdx rowIdx row index array for specifying mother parity check matrix.
         * @param rows_to_combine_rate_adapt rows_to_combine_rate_adapt array of marix line indices to be combined for rate adaption
         * @param initial_row_combs initial_row_combs number of line indices to combine initially
         */
        RateAdaptiveCode(const std::vector<colptr_t> &colptr,
                         const std::vector<idx_t> &rowIdx,
                         const std::vector<idx_t> &rows_to_combine_rate_adapt,
                         std::size_t initial_row_combs = 0)
                : colptr(colptr), row_idx(rowIdx), rows_to_combine(rows_to_combine_rate_adapt),
                  n_mother_rows(*std::max_element(rowIdx.begin(), rowIdx.end()) + 1u),
                  n_cols(colptr.size() - 1),
                  n_ra_rows(*std::max_element(rowIdx.begin(), rowIdx.end()) + 1u - initial_row_combs / 2) {
            if (rows_to_combine.size() % 2 != 0) {
                throw std::domain_error("The number of rows to combine for rate adaption "
                                        "(size of argument array) is an odd number (expected even).");
            }

            if (initial_row_combs >= rows_to_combine.size() / 2) {
                throw std::domain_error("The number of desired initial row combinations for rate adaption "
                                        "is larger than the given array of lines to combine.");
            }

            recompute_pos_vn_cn(initial_row_combs);
//            if (do_elimination_check && has_var_node_eliminations()) {
//                throw std::domain_error("Given rate adaption implies variable node eliminations. "
//                                        "Rate adaption with eliminations degrades performance. Do not use!");
//            }
        }


        // ---------------------------------------------------------------------------------------------- public methods
        constexpr void encode_no_ra(const std::vector<Bit> &in, std::vector<Bit> &out) const {
            if (in.size() != n_cols) {
                DEBUG_MESSAGE("Encoder (encode_no_ra) received invalid input length.");  // TODO maybe use exception?
                return;
            }
            out.assign(n_mother_rows, 0);

            for (std::size_t col = 0; col < in.size(); col++)
                for (std::size_t j = colptr[col]; j < colptr[col + 1]; j++)
                    out[row_idx[j]] = xor_as_bools(out[row_idx[j]], in[col]);
        }

        /// does not change internal rate adaption state!
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

            // int8_t because value -1 is used to mark bits included into final output
            std::vector<std::int8_t> non_ra_encoding(n_mother_rows);

            // first encode without rate adaption
            for (std::size_t col = 0; col < in.size(); col++)
                for (std::size_t j = colptr[col]; j < colptr[col + 1]; j++)
                    non_ra_encoding[row_idx[j]] = xor_as_bools(non_ra_encoding[row_idx[j]], in[col]);

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

        /// compute log-likelihood-ratios for a given keye and channel parameter
        static std::vector<double> llrs_bsc(const std::vector<Bit> &bitstring, const double bsc_channel_parameter) {
            double vlog = log((1 - bsc_channel_parameter) / bsc_channel_parameter);
            std::vector<double> llrs(bitstring.size());
            for (std::size_t i{}; i < llrs.size(); ++i) {
                llrs[i] = vlog * (1 - 2 * bitstring[i]); // log likelihood ratios
            }
            return llrs;
        }

        /// decoder infers rate from the length of the syndrome and changes the internal decoder state to match this rate.
        /// Note: since this function modifies the code (by performing rate adaption), it is NOT CONST.
        /// this change may be somewhat computationally expensive TODO benchmark this
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
         * Belief propagation decoder
         *
         * @param llrs: Log likelihood ratios representing the received message
         * @param syndrome: Syndrome of the sent message
         * @param out: Buffer to which the function writes its prediction for the sent message.
         * @param max_num_iter: Maximum number of iterations for the PB algorithm.
         *      Note that the algorithm always terminates automatically when the current prediction matches
         *      the syndrome (early termination), which means that the actual number of iterations cannot be controlled.
         * @param vsat: Cut-off value for messages.
         * @return true if and only if the syndrome of buffer `out` matches given `syndrome`,
         *      i.e., in case the decoder converged.
         */
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
                for (const auto &m : msg_v) {
                    for (const auto &v : m) {
                        if (std::isnan(v)) {
                            // TODO maybe use exception?
                            DEBUG_MESSAGE("Decoder Diverged at iteration " << it_unused);
                            return false;
                        }
                    }
                }
            }

            return false;  // Decoding was not successful.
        }

        /// manually trigger rate adaption. In normal circumstances, the user does not need this function
        void set_rate(std::size_t n_line_combs) {
            recompute_pos_vn_cn(n_line_combs);
        }


        /// TODO make private
        // TODO benchmark and compare performance to using pos_cn, as well as csc format
        constexpr void encode_at_current_rate(
                const std::vector<Bit> &in, std::vector<Bit> &out) const {
            if (in.size() != n_cols) {
                DEBUG_MESSAGE("Encoder received invalid input length.");  // TODO maybe use exception?
                return;
            }

            out.assign(pos_varn.size(), 0);

            for (std::size_t i{}; i < pos_varn.size(); ++i) {
                for (auto &var_node : pos_varn[i]) {
                    out[i] = xor_as_bools(out[i], in[var_node]);
                }
            }
        }

        bool operator==(const RateAdaptiveCode &rhs) const {
            return colptr == rhs.colptr &&
                   row_idx == rhs.row_idx &&
                   rows_to_combine == rhs.rows_to_combine &&
                   pos_checkn == rhs.pos_checkn &&
                   pos_varn == rhs.pos_varn &&
                   n_mother_rows == rhs.n_mother_rows &&
                   n_cols == rhs.n_cols &&
                   n_ra_rows == rhs.n_ra_rows;
        }

        bool operator!=(const RateAdaptiveCode &rhs) const {
            return !(rhs == *this); // Note: do not "simplify". The call to == is required.
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
        [[nodiscard]] std::size_t get_n_rows_mother_matrix() const {
            return n_mother_rows;
        }

        /// Includes rate adaption. Access to internal state!
        [[nodiscard]] std::size_t get_n_rows_after_rate_adaption() const {
            return n_ra_rows;
        }

        [[nodiscard]] std::size_t getNCols() const {
            return n_cols;
        }

        [[nodiscard]] std::size_t get_max_ra_steps() const {
            return rows_to_combine.size() / 2;
        }

    private:   // -------------------------------------------------------------------------------------- private members

        constexpr static bool xor_as_bools(Bit lhs, Bit rhs) {
            return (static_cast<bool>(lhs) != static_cast<bool>(rhs));
        }

        void check_node_update(std::vector<std::vector<double>> &msg_c,
                               const std::vector<std::vector<double>> &msg_v,
                               const std::vector<Bit> &syndrome) const {
            double msg_part{};
            std::vector<std::size_t> mc_position(n_cols);

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
                        DEBUG_MESSAGE("Decoder went into the 'untested bit'!!");
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
            std::vector<std::size_t> mv_position(n_cols);

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
         * Recompute inner representation of rate adapted LDPC code (pos_varn and pos_cn).
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

                std::vector<std::vector<idx_t>> pos_varn_nora{};

                if (n_line_combs == 0) {
                    for (idx_t col = 0; col < n_cols; col++) {
                        for (std::size_t j = colptr[col]; j < colptr[col + 1u]; j++) {
                            pos_varn[row_idx[j]].push_back(col);
                        }
                    }
                } else {
                    // first compute non-rate adapted variable nodes and store in temporary vector
                    pos_varn_nora.resize(n_mother_rows);

                    for (idx_t col = 0; col < n_cols; col++) {
                        for (std::size_t j = colptr[col]; j < colptr[col + 1u]; j++) {
                            pos_varn_nora[row_idx[j]].push_back(col);
                        }
                    }

                    // put results of combined lines at the back of the new LDPC code
                    const std::size_t start_of_ra_part = n_mother_rows - 2 * n_line_combs;

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
        // colptr and row_idx define the mother matrix, which each rate adaption starts from.
        const std::vector<colptr_t> colptr;
        const std::vector<idx_t> row_idx;
        const std::vector<idx_t> rows_to_combine;  // stores specification of rate adaption

        // these two (pos_checkn and pos_varn) store the current rate adapted code.
        std::vector<std::vector<idx_t>> pos_checkn;  // Input check nodes to each variable node
        std::vector<std::vector<idx_t>> pos_varn;  // Input variable nodes to each check node

        const std::size_t n_mother_rows;
        const std::size_t n_cols;
        std::size_t n_ra_rows;

    };

}

#endif //LDPC4QKD_LDPC_MATRIX_HPP
