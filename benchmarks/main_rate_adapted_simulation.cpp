//
// Created by alice on 07.09.21.
//


// Standard library
#include <iostream>
#include <random>
#include <chrono>

// Project scope
#include "rate_adaptive_code.hpp"

// Automatically generated C++ code that contains the LDPC matrix.
#include "autogen_ldpc_matrix_csc.hpp"

namespace {
    double h2(double p) {
        return -p * ::log(p) - (1 - p) * log(1 - p);
    }


    template<typename T>
    void noise_bitstring_inplace(std::mt19937_64 &rng, std::vector<T> &src, double err_prob) {
        std::bernoulli_distribution distribution(err_prob);

        for (std::size_t i = 0; i < src.size(); i++) {
            if (distribution(rng)) {
                src[i] = !src[i];
            } else {
                src[i] = src[i];
            }
        }
    }


    LDPC4QKD::RateAdaptiveCode<bool> get_code_big_ra() {
        std::vector<std::uint32_t> colptr(AutogenLDPC::colptr.begin(), AutogenLDPC::colptr.end());
        std::vector<std::uint16_t> row_idx(AutogenLDPC::row_idx.begin(), AutogenLDPC::row_idx.end());
        return LDPC4QKD::RateAdaptiveCode<bool>(colptr, row_idx);
    }

}

void print_command_line_help() {
    std::cout << "Expecting exactly 5 arguments." << std::endl;
    std::cout << "Example arguments: <executable> 0.05 5000 100 50 42 200" << std::endl;
    std::cout << "Specifying:\n"
                 "BSC channel parameter\n"
                 "max. nr. of frames to test\n"
                 "nr. of frame errors at which to quit\n"
                 "max. number of BP algorithm iterations\n"
                 "Mersenne Twister seed\n"
                 "Update console output every n frames" << std::endl;
}


int main(int argc, char *argv[]) {
    std::mt19937_64 rng(42);
    auto H = get_code_big_ra();

    constexpr double p = 0.04;
    constexpr std::size_t num_frames_to_test = 5;
    constexpr std::uint8_t max_num_iter = 50;
    constexpr std::size_t rate_step = 10;

    std::vector<std::size_t> syndrome_size_success{};
    std::size_t frame_idx{1};  // counts the number of iterations
    for (; frame_idx < num_frames_to_test; ++frame_idx) {
        std::size_t current_syndrome_size = H.get_n_rows_mother_matrix();
        std::size_t success_syndrome_size = H.getNCols(); // assume whole codeword leaked unless decoding success

        for (; current_syndrome_size > H.get_max_ra_steps(); current_syndrome_size -= rate_step) {
            std::vector<bool> x(H.getNCols()); // true data sent over a noisy channel
            noise_bitstring_inplace(rng, x, 0.5);  // choose it randomly.

            std::vector<bool> syndrome;  // syndrome for error correction, which is sent over a noise-less channel.
            H.encode_with_ra(x, syndrome, current_syndrome_size);

            std::vector<bool> x_noised = x; // copy for distorted data
            noise_bitstring_inplace(rng, x_noised, p);

// log likelihood ratio (llr) computation
            double vlog = ::log((1 - p) / p);
            std::vector<double> llrs(x.size());
            for (std::size_t i{}; i < llrs.size(); ++i) {
                llrs[i] = vlog * (1 - 2 * x_noised[i]); // log likelihood ratios
            }

            std::vector<bool> solution;
            bool success = H.decode_infer_rate(llrs, syndrome, solution, max_num_iter);

            if (solution == x) {
                if (success) {
                    success_syndrome_size = syndrome.size();
                } else {
                    std::cerr << "DECODER GIVES CORRECT RESULT ALTHOUGH IT HAS NOT CONVERGED!!!!" << std::endl;
                }
            } else {
                if (success) {
                    std::cerr << "\n\nDECODER CONVERGED TO WRONG CODEWORD!!!!\n" << std::endl;
                }
            }
            syndrome_size_success.push_back(success_syndrome_size);
        }
    }

    std::cout << "all syndrome sizes:" << std::endl;
    for (auto s : syndrome_size_success) {
        std::cout << s << ' ';
    }
    std::cout << "\n\n";

    double avg_synd_size =
            static_cast<double>(std::accumulate(syndrome_size_success.begin(), syndrome_size_success.end(), 0ul))
            / static_cast<double>(syndrome_size_success.size());
    std::cout << "Average syndrome size (out of " << num_frames_to_test << " ): " << avg_synd_size << std::endl;

    double avg_rate = avg_synd_size / static_cast<double>(H.getNCols());
    std::cout << "Average rate: " << avg_rate << " (1.4 * h2(p) = " << 1.4 * h2(p) << ")" << std::endl;

}
