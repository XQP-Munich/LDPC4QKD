//
// Created by alice on 07.09.21.
// In a Slepian-Wolf coding setting, for a given codeword and noised codeword, there is a minimum coding rate at which
// the syndrome decoding succeeds. This program determines the average minimum coding rate across many noised codewords.
//


// Standard library
#include <iostream>
#include <random>
#include <chrono>

// Command line argument parser library
#include "CmdParser-1.1.0/cmdparser.hpp"

// Project scope
#include "rate_adaptive_code.hpp"

#include "code_simulation_helpers.hpp"

using namespace LDPC4QKD::CodeSimulationHelpers;


template<typename RateAdaptiveCodeTemplate>
std::vector<std::size_t> run_simulation(RateAdaptiveCodeTemplate &H,
                                        double p,
                                        std::size_t num_frames_to_test,
                                        std::mt19937_64 &rng,
                                        std::size_t max_num_iter = 50,
                                        std::size_t update_console_every_n_frames = 100) {
    // assume whole codeword leaked unless decoding success
    std::vector<std::size_t> succesful_syndrome_sizes(num_frames_to_test, H.getNCols());
    constexpr int ra_step_accuracy = 1;

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


void configure_parser(cli::Parser &parser) {
    parser.set_optional<std::size_t>(
            "s", "seed", 42,
            "Mersenne Twister seed. Used to generate random bit-strings and simulate the noise channel.");

    parser.set_optional<std::size_t>(
            "upn", "update-console-n-frames", 100,
            "Update console output every n frames");

    parser.set_optional<std::size_t>(
            "nf", "num-frames-to-test", 1,
            "Number of frames to test (find optimal rate for).");

    parser.set_optional<std::size_t>(
            "i", "iter-bp", 50,
            "Maximum number of belief propagation (BP) algorithm iterations.");

    parser.set_optional<std::size_t>(
            "me", "max-frame-errors", 50,
            "Number of frame errors at which to quit the simulation. Specify zero for 'no condition'.");

    parser.set_optional<double>(
            "p", "channel-parameter", 0.02,
            "Binary Symmetric Channel (BSC) channel parameter. I.e., probability of a bit to be flipped.");

    parser.set_required<std::string>(
            "cp", "code-path",
            "Path to file containing LDPC code (`.cscmat` format. Note: does not accept QC exponents!)");

    parser.set_required<std::string>(
            "rp", "rate-adaption-path",
            "Path to file containing rate adaption for the LDPC code (`csv` format. Two columns of indices).");
}


int main(int argc, char *argv[]) {
    // parse command line arguments
    cli::Parser parser(argc, argv);
    configure_parser(parser);
    parser.run_and_exit_if_error();

    auto p = parser.get<double>("p");
    auto num_frames_to_test = parser.get<std::size_t>("nf");
    auto max_bp_iter = parser.get<std::size_t>("i");
    auto rng_seed = parser.get<std::size_t>("s");
    auto update_console_every_n_frames = parser.get<std::size_t>("upn");
    auto cscmat_file_path = parser.get<std::string>("cp");;
    auto rate_adaption_file_path = parser.get<std::string>("rp");;

    auto H = load_ldpc(cscmat_file_path, rate_adaption_file_path);

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
            max_bp_iter, update_console_every_n_frames);

    auto now = std::chrono::steady_clock::now();
    std::cout << "\n\nDONE! Simulation time: " <<
              std::chrono::duration_cast<std::chrono::seconds>(now - begin).count() << " seconds." << '\n';


    std::cout << "all syndrome sizes:" << std::endl;
    for (auto s : syndrome_size_success) {
        std::cout << s << ' ';
    }
    std::cout << "\n\n";

    double avg_synd_size = avg(syndrome_size_success);
    std::cout << "Average syndrome size (out of " << num_frames_to_test << " codewords tried): " << avg_synd_size << std::endl;

    double avg_rate = avg_synd_size / static_cast<double>(H.getNCols());
    std::cout << "Average rate: " << avg_rate << " (inefficiency f = " << avg_rate / h2(p) << ")" << std::endl;
    exit(EXIT_SUCCESS);
}
