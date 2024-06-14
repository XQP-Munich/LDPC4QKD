//
// Created by alice on 23.04.21.
//

// This code file is a simple and self-contained (except for one header file that contains the decoder implementation)
// example showing how to use the LDPC belief-propagation decoder.

// Project scope
#include "LDPC4QKD/rate_adaptive_code.hpp"


/// Get a LDPC matrix
auto get_code_small() {
    /// We use this matrix as an example:
    ///    H =  [1 0 1 0 1 0 1
    ///			0 1 1 0 0 1 1
    ///			0 0 0 1 1 1 1]
    ///
    /// To use it, we must convert H to compressed sparse column (CSC) storage:
    std::vector<std::uint32_t> colptr{0, 1, 2, 4, 5, 7, 9, 12};
    std::vector<std::uint16_t> row_idx{0, 1, 0, 1, 2, 0, 2, 1, 2, 0, 1, 2};
    return LDPC4QKD::RateAdaptiveCode(colptr, row_idx);
}


int main() {
    auto H = get_code_small();  // Get a small toy-example LDPC matrix

    // if you don't like `vector<bool>`, you can use e.g. `vector<uint8>` and store bits using only values 0 or 1.
    std::vector<bool> x{1, 1, 1, 1, 0, 0, 0}; // true data to be sent

    // syndrome computation using the LDPC matrix (sparse matrix-vector product modulo 2).
    std::vector<bool> syndrome;
    H.encode_no_ra(x, syndrome);

    std::vector<bool> x_noised{1, 1, 1, 1, 0, 0, 1}; // distorted data
    double p = 1. / 7; // channel error probability (for x_noised, we flipped 1 symbol out of 7)

    // Initial message computation for the decoder.
    // The decoder expects log-likelihood ratios (llr's) as input.
    double vlog = log((1 - p) / p);
    std::vector<double> llrs(x.size());
    for (std::size_t i{}; i < llrs.size(); ++i) {
        llrs[i] = vlog * (1 - 2 * x_noised[i]); // log likelihood ratios
    }
    // alternatively, use the built-in convenience method that does the same:
//    std::vector<double> llrs = LDPC4QKD::RateAdaptiveCode::llrs_bsc(x, p);

    std::vector<bool> solution;
    bool decoder_converged = H.decode_at_current_rate(llrs, syndrome, solution);

    bool decoding_correct = (solution == x);

    // return 0 if decoding was successful and correct, otherwise return 1.
    return !(decoder_converged & decoding_correct);
}
