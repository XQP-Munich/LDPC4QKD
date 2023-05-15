//
// Created by alice on 09.06.21.
// Note: Names and meaning of command line parameters are defined below.
//
constexpr auto help_text =
        "Frame Error Rate (FER) Simulator for Rate Adapted LDPC Codes\n"
        "\n"
        "This software is used to \n"
        "- load an LDPC code (from a .cscmat or bincsc.json file storing the full binary LDPC matrix in compressed sparse column (CSC) format, no QC exponents allowed!)\n"
        "- load rate adaption (from a csv file, list of pairs of row indices combined at each rate adaption step) "
        "   (this is optional; without rate adaption, only FER of the LDPC code can be simulated)\n"
        "- Simulate the FER of the given LDPC code at specified amount of rate adaption.";

// Standard library
#include <iostream>
#include <random>
#include <chrono>

// Command line argument parser library
#include "external/CmdParser-91aaa61e/cmdparser.hpp"

// Project scope
#include "rate_adaptive_code.hpp"
#include "code_simulation_helpers.hpp"

using namespace LDPC4QKD::CodeSimulationHelpers;


template<typename colptr_t=std::uint32_t, // integer type that fits ("number of non-zero matrix entries" + 1)
        typename idx_t=std::uint16_t>
std::pair<size_t, size_t> run_simulation(
        const LDPC4QKD::RateAdaptiveCode<bool, colptr_t, idx_t> &H,
        double p,
        std::size_t num_frames_to_test,
        std::mt19937_64 &rng,
        std::size_t max_num_iter = 50,
        std::size_t update_console_every_n_frames = 100,
        std::size_t quit_at_n_errors = 100) {
    std::size_t num_frame_errors{};
    std::size_t it{1};  // counts the number of iterations
    for (; num_frames_to_test == 0 || it < num_frames_to_test + 1; ++it) {
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
        if (update_console_every_n_frames && it % update_console_every_n_frames == 0) {
            std::cout << "current: " << num_frame_errors << " frame errors out of " << it
                      << " (FER~" << static_cast<double>(num_frame_errors) / static_cast<double>(it)
                      << ")..." << std::endl;
        }
        if (quit_at_n_errors != 0 && num_frame_errors >= quit_at_n_errors) {
            std::cout << "Quit simulation as max number of frame errors was reached." << std::endl;
            return std::make_pair(num_frame_errors, it);
        }
    }

    // minus 1 because the loop increments one more than it simulates.
    return std::make_pair(num_frame_errors, it - 1);
}


void configure_parser(cli::Parser &parser) {
    parser.set_optional<std::size_t>(
            "s", "seed", 42,
            "Mersenne Twister seed. Used to generate random bit-strings and simulate the noise channel.");

    parser.set_optional<std::size_t>(
            "upn", "update-console-n-frames", 100,
            "Update console output every n frames");

    parser.set_optional<std::size_t>(
            "mf", "max-frames", 0,
            "Maximum number of frames to test. Other conditions may terminate the simulation.");

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
            "Path to file containing LDPC code (`.cscmat` or `bincsc.json` format. Note: does not accept QC exponents!)");

    parser.set_optional<std::string>(
            "rp", "rate-adaption-path", "",
            "Path to file containing rate adaption for the LDPC code (`csv` format. Two columns of indices). "
            "If unspecified, no rate adaption is available.");

    parser.set_optional<std::size_t>(
            "rn", "rate-adaption-steps", 0,
            "Amount of rate adaption (number of row combinations) used for the simulation."
            "Can only be non-zero if a rate adaption file is also given.");
}


int main(int argc, char *argv[]) {
    // parse command line arguments
    cli::Parser parser(argc, argv, help_text);
    configure_parser(parser);
    parser.run_and_exit_if_error();

    auto p = parser.get<double>("p");
    auto max_num_frames_to_test = parser.get<std::size_t>("mf");
    auto quit_at_n_errors = parser.get<std::size_t>("me");
    auto max_bp_iter = parser.get<std::size_t>("i");
    auto rng_seed = parser.get<std::size_t>("s");
    auto update_console_every_n_frames = parser.get<std::size_t>("upn");
    auto code_file_path = parser.get<std::string>("cp");
    auto rate_adaption_file_path = parser.get<std::string>("rp");
    auto n_line_combs = parser.get<std::size_t>("rn");

    // create LDPC code, with rate adaption if specified.
    auto H = load_ldpc(code_file_path, rate_adaption_file_path);
    // set rate adaption. Only works if rate adaption was specified!
    H.set_rate(n_line_combs);

    // print received arguments (simulation parameters)
    std::cout << std::endl;
    std::cout << "Code path: '" << code_file_path << "'\n";
    std::cout << "Rate adaption path: '" << rate_adaption_file_path << "'\n";
    std::cout << "Running FER decoding test on channel parameter p : " << p << '\n';
    std::cout << "Max number of BP decoder iterations: " << static_cast<int>(max_bp_iter) << '\n';
    std::cout << "Max number of frames to simulate: " << max_num_frames_to_test << '\n';
    std::cout << "Quit at n frame errors: " << quit_at_n_errors << '\n';
    std::cout << "PRNG seed: " << rng_seed << '\n';
    std::cout << "Update console every n frames: " << update_console_every_n_frames << '\n';
    std::cout << "Code size before rate adaption: " << H.get_n_rows_mother_matrix() << " x " << H.getNCols() << '\n';
    std::cout << "Code size after rate adaption (if applicable): "
              << H.get_n_rows_after_rate_adaption() << " x " << H.getNCols() << "\n\n" << std::endl;

    std::mt19937_64 rng(rng_seed);
    auto begin = std::chrono::steady_clock::now();

    // perform frame error rate simulation.
    std::pair<std::size_t, std::size_t> result = run_simulation(H, p, max_num_frames_to_test, rng,
                                                                max_bp_iter,
                                                                update_console_every_n_frames, quit_at_n_errors);
    std::size_t num_frame_errors = result.first;
    std::size_t num_frames_tested = result.second;
    double naive_fer = static_cast<double>(num_frame_errors) / static_cast<double>(num_frames_tested);

    auto now = std::chrono::steady_clock::now();
    std::cout << "\n\nDONE! Simulation time: " <<
              std::chrono::duration_cast<std::chrono::seconds>(now - begin).count() << " seconds." << '\n';

    std::cout << "Recorded " << num_frame_errors << " frame errors out of " << num_frames_tested
              << " (FER~" << naive_fer << ")..." << std::endl;

    exit(EXIT_SUCCESS);
}
