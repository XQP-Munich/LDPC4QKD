//
// Created by alice on 07.09.21.
//
constexpr auto help_text =
        "Slepian-Wolf Critical Code Rate Simulator for Rate Adapted LDPC Codes\n"
        "\n"
        "In a Slepian-Wolf coding setting, for a given codeword and noised codeword, there is a minimum coding rate at which the syndrome decoding succeeds.\n"
        "This program determines the average minimum coding rate across many noised codewords.";

// Standard library
#include <iostream>
#include <random>
#include <chrono>

// Command line argument parser library
#include "cmdparser.hpp"

// Project scope
#undef LDPC4QKD_DEBUG_MESSAGES_ENABLED

#include "LDPC4QKD/rate_adaptive_code.hpp"
#include "LDPC4QKD/encoder_advanced.hpp"

#include "code_simulation_helpers.hpp"

using namespace LDPC4QKD::CodeSimulationHelpers;

template<typename T>
void print_avg(std::span<T> vals) {
    auto sum = std::accumulate(vals.begin(),
                               vals.end(), static_cast<T>(0));
    std::cout << "Simulated " << vals.size() << " frames. Current average final syndrome size: "
              << static_cast<double>(sum) / static_cast<double>(vals.size()) << std::endl;
}


template<typename RateAdaptiveCodeTemplate>
std::vector<std::size_t> run_simulation_nobisect(RateAdaptiveCodeTemplate &H,
                                                 const double p,
                                                 const std::size_t num_frames_to_test,
                                                 std::mt19937_64 &rng,
                                                 const std::size_t step_size,
                                                 const std::size_t start_syndrome_size,
                                                 const std::size_t max_num_iter = 50,
                                                 const std::size_t update_console_every_n_frames = 100) {
    // assume whole codeword leaked unless decoding success
    std::vector<std::size_t> final_syndrome_sizes(num_frames_to_test, H.getNCols());
    std::vector<std::size_t> n_messages(num_frames_to_test, 0);

    std::size_t frame_idx{0};  // counts the number of iterations
    for (; frame_idx < num_frames_to_test; ++frame_idx) {

        std::size_t current_syndrome_size = start_syndrome_size;

        while (true) {
            std::vector<bool> x(H.getNCols()); // true data sent over a noisy channel
            noise_bitstring_inplace(rng, x, 0.5);  // choose it randomly.

            std::vector<bool> syndrome;  // syndrome for error correction, which is sent over a noise-less channel.
            H.encode_with_ra(x, syndrome, current_syndrome_size);
            ++(n_messages.at(frame_idx));

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
                final_syndrome_sizes.at(frame_idx) = syndrome.size();
                break;  // Go to next codeword!
            } else {
                if (success) {
                    // here we stop adding new syndrome because decoder thinks it succeeded.
                    // The error is noticed later, but we failed overall and record worst-case syndrome size.
                    // This assignment is redundant (default value is written again) but we want to be explicit.
                    final_syndrome_sizes.at(frame_idx) = H.getNCols();
                    break;
                }
                current_syndrome_size += step_size; //  try again with more syndrome.

            }
        }
        if (update_console_every_n_frames && frame_idx % update_console_every_n_frames == 0) {
            print_avg(std::span{final_syndrome_sizes.begin(),
                                std::next(final_syndrome_sizes.begin(), static_cast<int>(frame_idx + 1))});
            for (auto x: std::span{n_messages.begin(),
                                   std::next(n_messages.begin(), static_cast<int32_t>(frame_idx + 1))}) {
                std::cout << " " << x;
            }
            std:: cout << " ." << '\n';
        }
    }
    std::cout << std::endl;

    return final_syndrome_sizes;
}


template<typename RateAdaptiveCodeTemplate>
std::vector<std::size_t> run_simulation_bisect(RateAdaptiveCodeTemplate &H,
                                               const double p,
                                               const std::size_t num_frames_to_test,
                                               std::mt19937_64 &rng,
                                               const std::size_t max_num_iter = 50,
                                               const std::size_t update_console_every_n_frames = 100) {
    // assume whole codeword leaked unless decoding success
    std::vector<std::size_t> final_syndrome_sizes(num_frames_to_test, H.getNCols());

    constexpr int ra_step_accuracy = 1;

    std::size_t frame_idx{0};  // counts the number of iterations
    for (; frame_idx < num_frames_to_test; ++frame_idx) {

        std::size_t min_syndrome_size = H.get_n_rows_mother_matrix() - H.get_max_ra_steps();
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
                final_syndrome_sizes.at(frame_idx) = syndrome.size();
                max_syndrome_size = syndrome.size();
            } else {
                if (success) {
                    std::cerr << "DECODER CONVERGED TO WRONG CODEWORD!!!!" << std::endl;
                }
                min_syndrome_size = syndrome.size();
            }
        }
        if (update_console_every_n_frames && frame_idx % update_console_every_n_frames == 0) {
            print_avg(std::span{final_syndrome_sizes.begin(),
                                std::next(final_syndrome_sizes.begin(), static_cast<int>(frame_idx + 1))});
        }
    }
    std::cout << std::endl;

    return final_syndrome_sizes;
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

    parser.set_optional<std::size_t>(
            "d", "step-size", 0,
            "If non-zero, we start at start-syndrome-size and increase the syndrome in steps of step-size.");

    parser.set_optional<std::size_t>(
            "s", "start-syndrome-size", 0,
            "Only has an effect if step-size is non-zero. "
            "If it has an effect, start at this syndrome size and increase until success.");

    parser.set_optional<double>(
            "p", "channel-parameter", 0.02,
            "Binary Symmetric Channel (BSC) channel parameter. I.e., probability of a bit to be flipped.");

    parser.set_required<std::string>(
            "cp", "code-path",
            "Path to file containing LDPC code (`.cscmat` or `.bincsc.json` format. Note: does not accept QC exponents!)."
            "Alternatively, a numeric index selecting pre-compiled code (e.g., `3` selects code 2048x4096 hash 0c809c3).");

    parser.set_optional<std::string>(
            "rp", "rate-adaption-path", "",
            "Path to file containing rate adaption for the LDPC code (`csv` format. Two columns of indices)."
            "This is ignored (and should be left to default) if `code-path` is specified as a numeric index.");
}


int main(int argc, char *argv[]) {
    // parse command line arguments
    cli::Parser parser(argc, argv, help_text);
    configure_parser(parser);
    parser.run_and_exit_if_error();

    const auto p = parser.get<double>("p");
    const auto num_frames_to_test = parser.get<std::size_t>("nf");
    const auto max_bp_iter = parser.get<std::size_t>("i");
    const auto rng_seed = parser.get<std::size_t>("s");
    const auto update_console_every_n_frames = parser.get<std::size_t>("upn");
    const auto code_file_path = parser.get<std::string>("cp");;
    const auto rate_adaption_file_path = parser.get<std::string>("rp");

    std::cout << std::endl;
    // try to parse code path as index. If not possible, it's actually a path, and then we load from that file.
    LDPC4QKD::RateAdaptiveCode<std::uint32_t> H = [&code_file_path, &rate_adaption_file_path]() {
        try {
            int idx = std::stoi(code_file_path); // we expect this to fail.
            auto H = LDPC4QKD::HelperFixedSize::get_rate_adaptive_code(static_cast<size_t>(idx));
            std::cout << "Pre-compiled LDPC Code and rate adaption used from index: " << idx << '\n';
            if (!rate_adaption_file_path.empty()) {
                std::cerr << "Warning: given rate adaption file path is ignored!" << std::endl;
            }
            return H;
        } catch (...) {
            auto H = load_ldpc(code_file_path, rate_adaption_file_path);
            std::cout << "LDPC Code loaded from file: " << code_file_path << '\n';
            std::cout << "Rate adaption loaded from file: " << rate_adaption_file_path << '\n';
            return H;
        }
    }();

    const auto step_size = parser.get<std::size_t>("d");
    const bool run_bisection = (step_size == 0);
    const auto min_supported_syndrome_size = H.get_n_rows_mother_matrix() - H.get_max_ra_steps();
    const auto start_syndrome_size = [&parser, &min_supported_syndrome_size]() {
        auto user_requested_start_syndrome_size = parser.get<std::size_t>("s");
        if (user_requested_start_syndrome_size < min_supported_syndrome_size) {
            std::cerr << "Warning:  Requested start syndrome size not supported. "
                         "Using smallest supported by selected code, which is "
                      << min_supported_syndrome_size << std::endl;
            return min_supported_syndrome_size;
        }
        return user_requested_start_syndrome_size;
    }();

    if (run_bisection && start_syndrome_size != 0) {
        std::cerr << "Warning: ignoring start-syndrome-size because bisection rate search is enabled "
                     "(where we check all available rates anyway)."
                  << std::endl;
    }

    std::cout << "Code size: " << H.get_n_rows_after_rate_adaption() << " x " << H.getNCols() << '\n';
    std::cout << "Running FER decoding test on channel parameter p : " << p << '\n';
    std::cout << "Max number decoder iterations: " << static_cast<int>(max_bp_iter) << '\n';
    std::cout << "Number of frames to simulate: " << num_frames_to_test << '\n';
    std::cout << "PRNG seed: " << rng_seed << '\n';
    if (run_bisection) {
        std::cout << "Running bisection on syndrome sizes from " << H.get_n_rows_mother_matrix() - H.get_max_ra_steps()
                  << " to " << H.get_n_rows_mother_matrix() << '\n';
    } else {
        std::cout << "step-size: " << step_size << '\n';
        std::cout << "start-syndrome-size: " << start_syndrome_size << '\n';
    }
    std::cout << "\n" << std::endl;

    std::mt19937_64 rng(rng_seed);
    auto begin = std::chrono::steady_clock::now();

    const auto syndrome_size_success = [&]() {
        if (run_bisection) {
            return run_simulation_bisect(
                    H, p, num_frames_to_test, rng,
                    max_bp_iter, update_console_every_n_frames);
        } else {
            return run_simulation_nobisect(
                    H, p, num_frames_to_test, rng, step_size, start_syndrome_size,
                    max_bp_iter, update_console_every_n_frames);
        }
    }();

    auto now = std::chrono::steady_clock::now();
    std::cout << "\n\nDONE! Simulation time: " <<
              std::chrono::duration_cast<std::chrono::seconds>(now - begin).count() << " seconds." << '\n';


    std::cout << "all syndrome sizes:" << std::endl;
    for (auto s: syndrome_size_success) {
        std::cout << s << ' ';
    }
    std::cout << "\n\n";

    double avg_synd_size = avg(syndrome_size_success);
    std::cout << "Average syndrome size (out of " << num_frames_to_test << " codewords tried): " << avg_synd_size
              << std::endl;

    double avg_rate = avg_synd_size / static_cast<double>(H.getNCols());
    std::cout << "Average rate: " << avg_rate << " (inefficiency f = " << avg_rate / h2(p) << ")" << std::endl;
    exit(EXIT_SUCCESS);
}
