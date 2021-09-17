//
// Created by alice on 09.06.21.
//

// Standard library
#include <iostream>
#include <random>
#include <chrono>

// Project scope
#include "rate_adaptive_code.hpp"

// Automatically generated C++ code that contains the LDPC matrix.
#include "autogen_ldpc_matrix_csc.hpp"


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


LDPC4QKD::RateAdaptiveCode<bool> get_code_big_nora() {
    std::vector<std::uint32_t> colptr(AutogenLDPC::colptr.begin(), AutogenLDPC::colptr.end());
    std::vector<std::uint16_t> row_idx(AutogenLDPC::row_idx.begin(), AutogenLDPC::row_idx.end());
    return LDPC4QKD::RateAdaptiveCode<bool>(colptr, row_idx);
}


std::pair<size_t, size_t> run_simulation(
        const LDPC4QKD::RateAdaptiveCode<bool> &H,
        double p,
        std::size_t num_frames_to_test,
        std::mt19937_64 &rng,
        std::uint8_t max_num_iter = 50,
        long update_console_every_n_frames = 100,
        long quit_at_n_errors = 100) {
    std::size_t num_frame_errors{};
    std::size_t it{1};  // counts the number of iterations
    for (; it < num_frames_to_test; ++it) {
        std::vector<bool> x(H.getNCols()); // true data sent over a noisy channel
        noise_bitstring_inplace(rng, x, 0.5);  // choose it randomly.

        std::vector<bool> syndrome;  // syndrome for error correction, which is sent over a noise-less channel.
        H.encode_no_ra(x, syndrome);

        std::vector<bool> x_noised = x; // distorted data
        noise_bitstring_inplace(rng, x_noised, p);

        double vlog = log((1 - p) / p);
        std::vector<double> llrs(x.size());
        for (std::size_t i{}; i < llrs.size(); ++i) {
            llrs[i] = vlog * (1 - 2 * x_noised[i]); // log likelihood ratios
        }

        std::vector<bool> solution;
        bool success = H.decode_at_current_rate(llrs, syndrome, solution, max_num_iter);

        if (success) {
            if (solution != x)
                std::cerr << "\n\nDECODER CONVERGED TO WRONG CODEWORD!!!!\n" << std::endl;
        } else {
            num_frame_errors++;
            if (solution == x)
                std::cerr << "DECODER GIVES CORRECT RESULT ALTHOUGH IT HAS NOT CONVERGED!!!!" << std::endl;
        }
        if (it % update_console_every_n_frames == 0) {
            std::cout << "current: " << num_frame_errors << " frame errors out of " << it
                      << " (FER~" << static_cast<double>(num_frame_errors) / it << ")..." << std::endl;
        }
        if (num_frame_errors >= quit_at_n_errors) {
            std::cout << "Quit simulation as max number of frame errors was reached." << std::endl;
            return std::make_pair(num_frame_errors, it);
        }
    }

    return std::make_pair(num_frame_errors, it);
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
    std::cout << "Program call:" << argv[0] << std::endl;

    std::vector<std::string> args{argv + 1, argv + argc};
    if (args.size() != 6) {
        std::cout << "Received " << args.size() << " arguments.\n" << std::endl;
        print_command_line_help();
        exit(EXIT_FAILURE);
    }

    double p{};
    std::size_t max_num_frames_to_test{};
    long quit_at_n_errors{};
    std::uint8_t max_bp_iter{};
    std::size_t rng_seed{};
    long update_console_every_n_frames{};

    try {
        p = stod(args[0]); // channel error probability
        max_num_frames_to_test = stol(args[1]);
        quit_at_n_errors = stol(args[2]);
        max_bp_iter = stoi(args[3]);
        rng_seed = stol(args[4]);
        update_console_every_n_frames = stol(args[5]);
    }
    catch (...) {
        std::cout << "Invalid command line arguments." << std::endl;
        print_command_line_help();
        exit(EXIT_FAILURE);
    }

    auto H = get_code_big_nora();

    std::cout << "Code size: " << H.get_n_rows_after_rate_adaption() << " x " << H.getNCols() << '\n';
    std::cout << "Running FER decoding test on channel parameter p : " << p << '\n';
    std::cout << "Max number decoder iterations: " << static_cast<int>(max_bp_iter) << '\n';
    std::cout << "Number of frames to simulate: " << max_num_frames_to_test << '\n';
    std::cout << "Quit at n frame errors: " << quit_at_n_errors << '\n';
    std::cout << "PRNG seed: " << rng_seed << '\n';
    std::cout << "\n" << std::endl;

    std::mt19937_64 rng(rng_seed);

    auto begin = std::chrono::steady_clock::now();

    std::pair<std::size_t, std::size_t> result = run_simulation(H, p, max_num_frames_to_test, rng,
                                                                max_bp_iter, update_console_every_n_frames);
    std::size_t num_frame_errors = result.first;
    std::size_t num_frames_tested = result.second;

    auto now = std::chrono::steady_clock::now();

    std::cout << "\n\nDONE! Simulation time: " <<
              std::chrono::duration_cast<std::chrono::seconds>(now - begin).count() << " seconds." << '\n';
    std::cout << "Recorded " << num_frame_errors << " frame errors out of " << num_frames_tested
              << " (FER~" << static_cast<double>(num_frame_errors) / num_frames_tested << ")..." << std::endl;

    exit(EXIT_SUCCESS);
}