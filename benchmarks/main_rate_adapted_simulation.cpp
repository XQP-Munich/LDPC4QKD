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
#include "autogen_rate_adaption.hpp"

#include "code_simulation_helpers.hpp"

using namespace LDPC4QKD::CodeSimulationHelpers;


void print_command_line_help() {
    std::cout << "Expecting exactly 7 arguments." << std::endl;
    std::cout << "Example arguments: <executable> 0.05 10 50 42 2 ./ldpc_filename.cscmat ./rate_adaption_filename.csv"
              << std::endl;
    std::cout << "Specifying:\n"
                 "BSC channel parameter\n"
                 "nr. of frames to test\n"
                 "max. number of BP algorithm iterations\n"
                 "Mersenne Twister seed\n"
                 "Update console output every n frames\n"
                 "Path to cscmat file containing LDPC code (not QC exponents!)\n"
                 "Path to csv file defining the rate adaption." << std::endl;
}

template<typename RateAdaptiveCodeTemplate>
std::vector<std::size_t> run_simulation(RateAdaptiveCodeTemplate &H,
                                        double p,
                                        std::size_t num_frames_to_test,
                                        std::mt19937_64 &rng,
                                        std::uint16_t max_num_iter = 50,
                                        long update_console_every_n_frames = 100,
                                        const int ra_step_accuracy = 1) {
    // assume whole codeword leaked unless decoding success
    std::vector<std::size_t> succesful_syndrome_sizes(num_frames_to_test, H.getNCols());

    std::size_t frame_idx{0};  // counts the number of iterations
    for (; frame_idx < num_frames_to_test; ++frame_idx) {
        std::size_t success_syndrome_size = H.getNCols(); // assume whole codeword leaked unless decoding success

        std::size_t min_syndrome_size = H.get_n_rows_mother_matrix() - H.get_max_ra_steps();;
        std::size_t max_syndrome_size = H.get_n_rows_mother_matrix();
        // bisection
        while ((max_syndrome_size - min_syndrome_size) > ra_step_accuracy) {
                        std::vector<bool> x(H.getNCols()); // true data sent over a noisy channel
            noise_bitstring_inplace(rng, x, 0.5);  // choose it randomly.

            std::vector<bool> syndrome;  // syndrome for error correction, which is sent over a noise-less channel.
            const std::size_t current_syndrome_size = (max_syndrome_size + min_syndrome_size) / 2;
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

            if (success && solution == x) {
                succesful_syndrome_sizes.at(frame_idx) = syndrome.size();
                max_syndrome_size = syndrome.size();
            } else {
                if (success) {
                    std::cerr << "\n\nDECODER CONVERGED TO WRONG CODEWORD!!!!\n" << std::endl;
                }
                min_syndrome_size = syndrome.size();
            }
        }
        if (update_console_every_n_frames && frame_idx % update_console_every_n_frames == 0) {
            std::cout << "\rcurrent average successful syndrome size: " << avg(succesful_syndrome_sizes);
        }
    }
    std::cout << std::endl;

    return succesful_syndrome_sizes;
}


int main(int argc, char *argv[]) {
    std::cout << "Program call: " << argv[0] << std::endl;
    std::vector<std::string> args{argv + 1, argv + argc};
    if (args.size() != 7) {
        std::cout << "Received " << args.size() << " arguments.\n" << std::endl;
        print_command_line_help();
        exit(EXIT_FAILURE);
    }

    constexpr std::size_t rate_step = 10;
    double p{};
    std::size_t num_frames_to_test{};
    std::uint16_t max_bp_iter{};
    std::size_t rng_seed{};
    long update_console_every_n_frames{};
    std::string cscmat_file_path;
    std::string rate_adaption_file_path;

    try {
        p = stod(args[0]); // channel error probability
        num_frames_to_test = stol(args[1]);
        max_bp_iter = stoi(args[2]);
        rng_seed = stol(args[3]);
        update_console_every_n_frames = stol(args[4]);
        cscmat_file_path = args[5];
        rate_adaption_file_path = args[6];
    }
    catch (...) {
        std::cout << "Invalid command line arguments." << std::endl;
        print_command_line_help();
        exit(EXIT_FAILURE);
    }

    auto H = get_code_big_wra(cscmat_file_path, rate_adaption_file_path);

    std::cout << std::endl;
    std::cout << "LDPC Code loaded from file: " << cscmat_file_path << '\n';
    std::cout << "Rate adaption loaded from file: " << rate_adaption_file_path << '\n';
    std::cout << "Code size: " << H.get_n_rows_after_rate_adaption() << " x " << H.getNCols() << '\n';
    std::cout << "Running FER decoding test on channel parameter p : " << p << '\n';
    std::cout << "Max number decoder iterations: " << static_cast<int>(max_bp_iter) << '\n';
    std::cout << "Number of frames to simulate: " << num_frames_to_test << '\n';
    std::cout << "PRNG seed: " << rng_seed << '\n';
    std::cout << "\n" << std::endl;

    std::mt19937_64 rng(rng_seed);
    auto begin = std::chrono::steady_clock::now();

    auto syndrome_size_success = run_simulation(
            H, p, num_frames_to_test, rng,
            max_bp_iter, update_console_every_n_frames, rate_step);

    auto now = std::chrono::steady_clock::now();
    std::cout << "\n\nDONE! Simulation time: " <<
              std::chrono::duration_cast<std::chrono::seconds>(now - begin).count() << " seconds." << '\n';


    std::cout << "all syndrome sizes:" << std::endl;
    for (auto s : syndrome_size_success) {
        std::cout << s << ' ';
    }
    std::cout << "\n\n";

    double avg_synd_size = avg(syndrome_size_success);
    std::cout << "Average syndrome size (out of " << num_frames_to_test << " ): " << avg_synd_size << std::endl;

    double avg_rate = avg_synd_size / static_cast<double>(H.getNCols());
    std::cout << "Average rate: " << avg_rate << " (inefficiency f = " << avg_rate / h2(p) << ")" << std::endl;
    exit(EXIT_SUCCESS);
}
