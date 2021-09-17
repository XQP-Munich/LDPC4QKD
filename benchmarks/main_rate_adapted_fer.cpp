//
// Created by alice on 09.06.21.
//

// Standard library
#include <iostream>
#include <random>
#include <chrono>

// Project scope
#include "rate_adaptive_code.hpp"

#include "code_simulation_helpers.hpp"
using namespace LDPC4QKD::CodeSimulationHelpers;


void print_command_line_help() {
    std::cout << "Expecting exactly 9 arguments." << std::endl;
    std::cout << "Example arguments: <executable> 0.05 5000 100 50 42 200 ./filename.cscmat ./rate_adaption_filename.csv 1000" << std::endl;
    std::cout << "Specifying:\n"
                 "BSC channel parameter\n"
                 "max. nr. of frames to test\n"
                 "nr. of frame errors at which to quit\n"
                 "max. number of BP algorithm iterations\n"
                 "Mersenne Twister seed\n"
                 "Update console output every n frames\n"
                 "Path to cscmat file containing LDPC code (not QC exponents!)\n"
                 "Path to csv file defining the rate adaption.\n"
                 "Amount of rate adaption (number of row combinations)" << std::endl;
}


template <typename colptr_t=std::uint32_t, // integer type that fits ("number of non-zero matrix entries" + 1)
        typename idx_t=std::uint16_t>
std::pair<size_t, size_t> run_simulation(
        const LDPC4QKD::RateAdaptiveCode<bool, colptr_t, idx_t> &H,
        double p,
        std::size_t num_frames_to_test,
        std::mt19937_64 &rng,
        std::uint16_t max_num_iter = 50,
        long update_console_every_n_frames = 100,
        long quit_at_n_errors = 100) {
    std::size_t num_frame_errors{};
    std::size_t it{1};  // counts the number of iterations
    for (; it < num_frames_to_test + 1; ++it) {
        std::vector<bool> x(H.getNCols()); // true data sent over a noisy channel
        noise_bitstring_inplace(rng, x, 0.5);  // choose it randomly.

        std::vector<bool> syndrome;  // syndrome for error correction, which is sent over a noise-less channel.
        H.encode_at_current_rate(x, syndrome);

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
            if (solution != x) {
                std::cerr << "\n\nDECODER CONVERGED TO WRONG CODEWORD!!!!\n" << std::endl;
                num_frame_errors++;
            }
        } else {
            num_frame_errors++;
            if (solution == x)
                std::cerr << "DECODER GIVES CORRECT RESULT ALTHOUGH IT HAS NOT CONVERGED!!!!" << std::endl;
        }
        if (it % update_console_every_n_frames == 0) {
            std::cout << "current: " << num_frame_errors << " frame errors out of " << it
                      << " (FER~" << static_cast<double>(num_frame_errors) / static_cast<double>(it)
                      << ")..." << std::endl;
        }
        if (num_frame_errors >= quit_at_n_errors) {
            std::cout << "Quit simulation as max number of frame errors was reached." << std::endl;
            return std::make_pair(num_frame_errors, it);
        }
    }

    // minus 1 because the loop increments one more than it simulates.
    return std::make_pair(num_frame_errors, it - 1);
}


int main(int argc, char *argv[]) {
    std::cout << "Program call:" << argv[0] << std::endl;

    std::vector<std::string> args{argv + 1, argv + argc};
    if (args.size() != 9) {
        std::cout << "Received " << args.size() << " arguments.\n" << std::endl;
        print_command_line_help();
        exit(EXIT_FAILURE);
    }

    double p{};
    std::size_t max_num_frames_to_test{};
    long quit_at_n_errors{};
    std::uint16_t max_bp_iter{};
    std::size_t rng_seed{};
    long update_console_every_n_frames{};
    std::string cscmat_file_path;
    std::string rate_adaption_file_path;
    std::size_t n_line_combs{};

    try {
        p = stod(args[0]); // channel error probability
        max_num_frames_to_test = stol(args[1]);
        quit_at_n_errors = stol(args[2]);
        max_bp_iter = stoi(args[3]);
        rng_seed = stol(args[4]);
        update_console_every_n_frames = stol(args[5]);
        cscmat_file_path = args[6];
        rate_adaption_file_path = args[7];
        n_line_combs = stol(args[8]);
    }
    catch (...) {
        std::cout << "Invalid command line arguments." << std::endl;
        print_command_line_help();
        exit(EXIT_FAILURE);
    }

    auto H = get_code_big_wra(cscmat_file_path, rate_adaption_file_path);
    H.set_rate(n_line_combs);

    std::cout << std::endl;
    std::cout << "Code path: " << cscmat_file_path << '\n';
    std::cout << "Rate adaption path: " << rate_adaption_file_path << '\n';
    std::cout << "Code size (before rate adaption): " << H.get_n_rows_mother_matrix() << " x " << H.getNCols() << '\n';
    std::cout << "Code size (after rate adaption): " << H.get_n_rows_after_rate_adaption() << " x " << H.getNCols() << '\n';
    std::cout << "Running FER decoding test on channel parameter p : " << p << '\n';
    std::cout << "Max number of BP decoder iterations: " << static_cast<int>(max_bp_iter) << '\n';
    std::cout << "Max number of frames to simulate: " << max_num_frames_to_test << '\n';
    std::cout << "Quit at n frame errors: " << quit_at_n_errors << '\n';
    std::cout << "PRNG seed: " << rng_seed << '\n';
    std::cout << "Update console every n frames: " << update_console_every_n_frames << '\n';
    std::cout << "\n" << std::endl;

    std::mt19937_64 rng(rng_seed);
    auto begin = std::chrono::steady_clock::now();

    std::pair<std::size_t, std::size_t> result = run_simulation(H, p, max_num_frames_to_test, rng,
                                                                max_bp_iter,
                                                                update_console_every_n_frames, quit_at_n_errors);
    std::size_t num_frame_errors = result.first;
    std::size_t num_frames_tested = result.second;

    auto now = std::chrono::steady_clock::now();

    std::cout << "\n\nDONE! Simulation time: " <<
              std::chrono::duration_cast<std::chrono::seconds>(now - begin).count() << " seconds." << '\n';
    std::cout << "Recorded " << num_frame_errors << " frame errors out of " << num_frames_tested
              << " (FER~" << static_cast<double>(num_frame_errors) / num_frames_tested << ")..." << std::endl;

    exit(EXIT_SUCCESS);
}
