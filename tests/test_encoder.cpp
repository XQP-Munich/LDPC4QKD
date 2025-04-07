//
// Created by alice on 23.04.21.
//

// Google Test framework
#include <gtest/gtest.h>
#include "helpers_for_testing.hpp"

// Standard library
#include <iostream>
#include <vector>
#include <chrono>

// Project scope
#include "fortest_autogen_ldpc_matrix_csc.hpp"
#include "fortest_autogen_rate_adaption.hpp"

#include "LDPC4QKD/encoder.hpp"


using namespace LDPC4QKD;
using namespace HelpersForTests;


TEST(encoder, Hx_basic_demo) {
    std::array<Bit, AutogenLDPC::N> input{};
    std::array<Bit, AutogenLDPC::M> output{};

    std::cout << "------------------------\n";
    std::cout << "Input nonzeros:\n";
    vec_to_arr<Bit, AutogenLDPC::N>(get_bitstring(AutogenLDPC::N), input);
    // noise_bitstring_inplace(input, 0.5, seed);

    print_nz_inds(input);

    std::cout << "------------------------\n";

    encode(input, output);

    std::cout << "Syndrome nonzeros:\n";
    print_nz_inds(output);

    std::cout << "------------------------\n";

    std::cout << "Syndrome hash: " << hash_vector(output) << std::endl;
}


TEST(encoder, Hx_benchmark) {
    std::array<Bit, AutogenLDPC::N> input{};
    std::array<Bit, AutogenLDPC::M> output{};

    constexpr int num_runs = 1000;
    unsigned int seed = 0;
    long long runtime = 0;

    std::cout << "Performing " << num_runs
        << " encodings for matrix size " << AutogenLDPC::M << " x " << AutogenLDPC::N << std::endl;
    std::cout << "Matrix has " << AutogenLDPC::num_nz << " non-zero entries." << std::endl;

    for (; seed < num_runs; ++seed) {
        noise_bitstring_inplace(input, 0.5, seed);

        auto begin = std::chrono::steady_clock::now();

        encode(input, output);

        auto now = std::chrono::steady_clock::now();

        auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(now - begin);
        runtime += elapsed.count();
        std::cout << "\rHash every run (makes sure the result is not optimized away): " << hash_vector(output);
    }
    std::cout << "\nAvg. runtime: " << static_cast<double>(runtime) / num_runs << " microseconds." << std::endl;
}


TEST(encoder, rate_adapted_Hx_benchmark) {
    constexpr std::size_t n_left_out_rate_adaption_steps = 2;
    std::array<Bit, AutogenLDPC::N> input{};
    std::array<Bit, AutogenLDPC::M> output{};
    std::array<Bit, AutogenLDPC::M/2 + n_left_out_rate_adaption_steps> ra_output{};

    constexpr int num_runs = 500;
    unsigned int seed = 0;
    long long runtime = 0;

    std::cout << "Performing " << num_runs
              << " encodings for matrix size " << AutogenLDPC::M << " x " << AutogenLDPC::N << std::endl;
    std::cout << "Matrix has " << AutogenLDPC::num_nz << " non-zero entries." << std::endl;

    for (; seed < num_runs; ++seed) {
        noise_bitstring_inplace(input, 0.5, seed);

        auto begin = std::chrono::steady_clock::now();

        encode(input, output);
        rate_adapt(output, ra_output);

        auto now = std::chrono::steady_clock::now();

        std::cout << "\rHash every run (makes sure the result is not optimized away): " << hash_vector(ra_output);

        auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(now - begin);
        runtime += elapsed.count();
    }
    std::cout << "\nAvg. runtime: " << static_cast<double>(runtime) / num_runs << " microseconds." << std::endl;
}


TEST(encoder, rate_adapted_Hx_benchmark_unsafe) {
    constexpr std::size_t n_left_out_rate_adaption_steps = 2;
    std::array<Bit, AutogenLDPC::N> input{};
    std::array<Bit, AutogenLDPC::M> output{};
    std::array<Bit, AutogenLDPC::M/2 + n_left_out_rate_adaption_steps> ra_output{};

    constexpr int num_runs = 500;
    unsigned int seed = 0;
    long long runtime = 0;

    std::cout << "Performing " << num_runs
              << " encodings for matrix size " << AutogenLDPC::M << " x " << AutogenLDPC::N << std::endl;
    std::cout << "Matrix has " << AutogenLDPC::num_nz << " non-zero entries." << std::endl;

    for (; seed < num_runs; ++seed) {
        noise_bitstring_inplace(input, 0.5, seed);

        auto begin = std::chrono::steady_clock::now();

        encode(input, output);
        rate_adapt(output, ra_output);

        auto now = std::chrono::steady_clock::now();

        std::cout << "\rHash every run (makes sure the result is not optimized away): " << hash_vector(ra_output);

        auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(now - begin);
        runtime += elapsed.count();
    }
    std::cout << "\nAvg. runtime: " << static_cast<double>(runtime) / num_runs << " microseconds." << std::endl;
}


TEST(encoder, rate_adapt_demo) {
    std::array<Bit, AutogenLDPC::N> input{
        1,1,1,1,0,0,0,0,1,1};
    std::array<Bit, AutogenLDPC::M> syndrome{};
    constexpr std::size_t rate_adaption_steps = AutogenRateAdapt::rows.size() / 2;
    std::array<Bit, AutogenLDPC::M - rate_adaption_steps> ratead_synd{};
//    vec_to_arr<Bit, AutogenLDPC::N>(get_bitstring(AutogenLDPC::N), input);

    std::cout << "------------------------\n";
    print_nz_inds(input);

    encode(input, syndrome);

    std::cout << "------------------------\n";
    print_nz_inds(syndrome);

    rate_adapt(syndrome, ratead_synd);

    std::cout << "------------------------\n";
    print_nz_inds(ratead_synd);
    std::cout << ratead_synd.size() << std::endl;
}


TEST(encoder, rate_adapt_unsafe_demo) {
    std::array<Bit, AutogenLDPC::N> input{
            1,1,1,1,0,0,0,0,1,1};
    std::array<Bit, AutogenLDPC::M> syndrome{};
    constexpr std::size_t rate_adaption_steps = AutogenRateAdapt::rows.size() / 2;
    std::array<Bit, AutogenLDPC::M - rate_adaption_steps> ratead_synd{};
//    vec_to_arr<Bit, AutogenLDPC::N>(get_bitstring(AutogenLDPC::N), input);

    std::cout << "------------------------\n";
    print_nz_inds(input);

    encode(input, syndrome);

    std::cout << "------------------------\n";
    print_nz_inds(syndrome);

    rate_adapt(syndrome, ratead_synd);

    auto result_copy = ratead_synd;
    ratead_synd = {};
    rate_adapt_unsafe(syndrome, ratead_synd.data(), ratead_synd.size());

    EXPECT_EQ(ratead_synd, result_copy);

    std::cout << "------------------------\n";
    print_nz_inds(ratead_synd);
    std::cout << ratead_synd.size() << std::endl;
}
