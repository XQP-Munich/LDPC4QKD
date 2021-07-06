#ifndef LDPC4QKD_LDPC_MATRIX_HPP
#define LDPC4QKD_LDPC_MATRIX_HPP

#include <cstdint>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>


#ifdef DEBUG_MESSAGES_ENABLED
#include <iostream>
#define DEBUG_MESSAGE(msg) do {std::cerr << msg << std::endl;} while (false)
#endif /* ifdef DEBUG_MESSAGES_ENABLED */


namespace LDPC4QKD {

    template<typename Bit, // type of matrix entires (only values zero and one are used, default to bool)
            typename colptr_t=std::uint32_t, // integer type that fits ("number of non-zero matrix entries" + 1)
            typename idx_t=std::uint16_t> // integer type fitting number of columns N (thus also number of rows M)
    class RateAdaptiveCode {
    public:
        // ------------------------------------------------------------------------------------------------ type aliases
        using NodeDegree = std::uint16_t;
        using MatrixIndex = idx_t;
        using ColumnPointer = colptr_t;
        using MatrixEntry = Bit;

        // ------------------------------------------------------------------------------------------------ constructors
        RateAdaptiveCode(const std::vector<colptr_t> &colptr, const std::vector<idx_t> &rowIdx)
                : colptr(colptr), row_idx(rowIdx),
                  n_cols(colptr.size() - 1),
                  n_rows(*std::max_element(rowIdx.begin(), rowIdx.end()) + 1) {
            recompute_check_node_degs();
            recompute_var_node_degs();
            recompute_pos_checkn();
            recompute_pos_varn();
        }

        // ---------------------------------------------------------------------------------------------- public methods
        constexpr void encode(const std::vector<Bit> &in, std::vector<Bit> &out) const {
            if (in.size() != n_cols) {
                DEBUG_MESSAGE("Encoder received invalid input length.");  // TODO maybe use exception?
                return;
            }
            out = std::vector<Bit>(n_rows);

            for (std::size_t col = 0; col < in.size(); col++)
                for (std::size_t j = colptr[col]; j < colptr[col + 1]; j++)
                    out[row_idx[j]] = xor_as_bools(out[row_idx[j]], in[col]);
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
        bool decode(const std::vector<double> &llrs,
                    const std::vector<Bit> &syndrome,
                    std::vector<Bit> &out,
                    const std::uint8_t max_num_iter = 50,
                    const double vsat = 100) const {
            // check inputs.
            if (llrs.size() != n_cols) {
                DEBUG_MESSAGE("Decoder received invalid input length."); // TODO maybe use exception?
                return false;
            }
            out.resize(llrs.size());

            std::vector<std::vector<double>> msg_v(n_rows);  // messages from variable nodes to check nodes
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
                encode(out, decision_syndrome);
                if (decision_syndrome == syndrome) {
                    return true;
                }

                // check for diverging decoder
                for (const auto &m : msg_v) {
                    for (const auto &v : m) {
                        if (std::isnan(v)) {
                            // TODO maybe use exception?
                            DEBUG_MESSAGE("Decoder Diverged at iteration" << it_unused);
                            return false;
                        }
                    }
                }
            }

            return false;  // Decoding was not successful.
        }

        // ----------------------------------------------------------------------------------------- getters and setters
        const std::vector<std::vector<idx_t>> &getPosCheckn() const {
            return pos_checkn;
        }

        const std::vector<std::vector<idx_t>> &getPosVarn() const {
            return pos_varn;
        }

        // TODO Unnecessary for users. Remove? Also: disallow any user access to colptr and row_idx.
        [[nodiscard]] colptr_t getNonzeros() const {
            return row_idx.size();
        }

        [[nodiscard]] std::size_t getNRows() const {
            return n_rows;
        }

        [[nodiscard]] std::size_t getNCols() const {
            return n_cols;
        }

        [[nodiscard]]
        const std::vector<NodeDegree> &getCheckNodeDegrees() const {
            return check_node_degrees;
        }

        [[nodiscard]]
        const std::vector<NodeDegree> &getVariableNodeDegrees() const {
            return variable_node_degrees;
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

            for (std::size_t m{}; m < n_rows; ++m) {
                // product of incoming messages
                double mc_prod = 1 - 2 * static_cast<double>(syndrome[m]);
                const auto curr_check_node_degree = check_node_degrees[m];
                for (std::size_t k{}; k < curr_check_node_degree; ++k) {
                    mc_prod *= ::tanh(0.5 * msg_v[m][k]);
                }

                for (std::size_t k{}; k < curr_check_node_degree; ++k) {
                    // computing message from
                    if (msg_v[m][k] == 0.) {  // TODO test this bit more carefully.
                        std::cerr << "Decoder went into the 'untested bit'!!" << std::endl;
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

                for (std::size_t k{}; k < variable_node_degrees[m]; ++k) {
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
        static void saturate(std::vector<std::vector<T>> mv, const double vsat) {
            for (auto &v : mv) {
                for (auto &a : v) {
                    if (a > vsat) { a = vsat; }
                    else if (a < -vsat) { a = -vsat; }
                }
            }
        }

        void recompute_check_node_degs() {
            check_node_degrees = std::vector<NodeDegree>(n_rows);

            for (std::size_t col = 0; col < n_cols; col++)
                for (std::size_t j = colptr[col]; j < colptr[col + 1]; j++)
                    check_node_degrees[row_idx[j]] += 1;
        }

        void recompute_var_node_degs() {
            variable_node_degrees = std::vector<NodeDegree>(n_cols);
            for (std::size_t col = 0; col < n_cols; col++)
                for (std::size_t j = colptr[col]; j < colptr[col + 1]; j++)
                    variable_node_degrees[col] += 1;
        }

        void recompute_pos_varn() {
            // This uses different size vectors for nodes with different degrees.
            // Alternatively, one could set the sizes to be the same using
//            auto max_check_deg = *std::max_element(check_node_degrees.begin(), check_node_degrees.end());
            pos_varn = std::vector<std::vector<idx_t>>(n_rows);

            for (std::size_t col = 0; col < n_cols; col++)
                for (std::size_t j = colptr[col]; j < colptr[col + 1]; j++)
                    pos_varn[row_idx[j]].push_back(col);
        }

        void recompute_pos_checkn() {
            // This uses different size vectors for nodes with different degrees.
            // Alternatively, one could set the sizes to be the same using
//            auto max_var_deg = *std::max_element(variable_node_degrees.begin(), variable_node_degrees.end());
            pos_checkn = std::vector<std::vector<idx_t>>(n_cols);

            for (std::size_t col = 0; col < n_cols; col++)
                for (std::size_t j = colptr[col]; j < colptr[col + 1]; j++)
                    pos_checkn[col].push_back(row_idx[j]);
        }

        // ---------------------------------------------------------------------------------------------- private fields
        std::vector<colptr_t> colptr;  // TODO declare const?
        std::vector<idx_t> row_idx;  // TODO declare const?

        std::vector<NodeDegree> check_node_degrees;  // TODO remove?
        std::vector<NodeDegree> variable_node_degrees;  // TODO remove?
        std::vector<std::vector<idx_t>> pos_checkn;  // Input check nodes to each variable node
        std::vector<std::vector<idx_t>> pos_varn;  // Input variable nodes to each check node

        std::size_t n_rows;
        std::size_t n_cols;
    };

}

#endif //LDPC4QKD_LDPC_MATRIX_HPP
