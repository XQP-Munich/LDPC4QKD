//
// Created by Adomas Baliuka on 02.05.24.
//

// Google Test framework
#include <gtest/gtest.h>
#include "helpers_for_testing.hpp"

// Standard library
#include <iostream>

// To be tested
#include "LDPC4QKD/prebuilt_codes.hpp"

using namespace LDPC4QKD;

namespace {
    void noise_bitstring_inplace(auto &src, double err_prob, unsigned int seed = 0) {
        std::mt19937_64 rng{seed}; // hard-coded seed for testing purposes.

        std::bernoulli_distribution distribution(err_prob);

        for (std::size_t i = 0; i < src.size(); i++) {
            if (distribution(rng)) {
                src[i] = !src[i];
            } else {
                src[i] = src[i];
            }
        }
    }
}

void debug_print_sizes() {
    std::cout
            // Protograph-based
            << "encoder_2048x6144_4663d91: " << sizeof(encoder_2048x6144_4663d91) << '\n'
            << "encoder_8192x24576_71b51c1: " << sizeof(encoder_8192x24576_71b51c1) << '\n'
            << "encoder_524288x1572864_4d78a9f: " << sizeof(encoder_524288x1572864_4d78a9f) << '\n'
            << "encoder_2048x4096_0c809c3: " << sizeof(encoder_2048x4096_0c809c3) << '\n'
            << "encoder_8192x16384_3fcad37: " << sizeof(encoder_8192x16384_3fcad37) << '\n'
            << "encoder_524288x1048576_9b50f98: " << sizeof(encoder_524288x1048576_9b50f98) << '\n'
            // degree-distribution-based
            << std::endl
            << "encoder_lrate_P1_block_819k: " << sizeof(encoder_lrate_P1_block_819k) << "\n"
            << "encoder_lrate_P15_block_819k: " << sizeof(encoder_lrate_P15_block_819k) << "\n"
            << "encoder_lrate_P2_block_819k: " << sizeof(encoder_lrate_P2_block_819k) << "\n"
            << "encoder_lrate_P25_block_819k: " << sizeof(encoder_lrate_P25_block_819k) << "\n"
            << "encoder_lrate_P3_block_819k: " << sizeof(encoder_lrate_P3_block_819k) << "\n"
            << "encoder_lrate_P35_block_819k: " << sizeof(encoder_lrate_P35_block_819k) << "\n"
            << "encoder_lrate_P4_block_819k: " << sizeof(encoder_lrate_P4_block_819k) << "\n"
            << "encoder_lrate_P45_block_819k: " << sizeof(encoder_lrate_P45_block_819k) << "\n"
            << "encoder_lrate_P5_block_819k: " << sizeof(encoder_lrate_P5_block_819k) << "\n";
}

TEST(test_prebuilt_codes, memory_usage_encoder_storage) {
    debug_print_sizes();

    constexpr std::size_t total_size_protograph =
            sizeof(encoder_2048x6144_4663d91) +
            sizeof(encoder_8192x24576_71b51c1) +
            sizeof(encoder_524288x1572864_4d78a9f) +
            sizeof(encoder_2048x4096_0c809c3) +
            sizeof(encoder_8192x16384_3fcad37) +
            sizeof(encoder_524288x1048576_9b50f98);

    constexpr std::size_t total_size_819k =
            sizeof(encoder_lrate_P1_block_819k) +
            sizeof(encoder_lrate_P15_block_819k) +
            sizeof(encoder_lrate_P2_block_819k) +
            sizeof(encoder_lrate_P25_block_819k) +
            sizeof(encoder_lrate_P3_block_819k) +
            sizeof(encoder_lrate_P35_block_819k) +
            sizeof(encoder_lrate_P4_block_819k) +
            sizeof(encoder_lrate_P45_block_819k) +
            sizeof(encoder_lrate_P5_block_819k);

    constexpr std::size_t total_size = total_size_protograph + total_size_819k;

    // Baselines measured directly (see this test's own `total_size`/`total_size_819k` computations)
    // at the time these assertions were added. Guards against silently ballooning embedded encoder
    // storage (e.g. from adding a new large prebuilt matrix, or a representation change) going
    // unnoticed. If a test fails because of a deliberate change, update the corresponding baseline
    // to the new measured value.
    constexpr std::size_t baseline_total_size = 258448;
    constexpr std::size_t baseline_819k_total_size = 156592;
    constexpr double margin = 1.01; // allow up to 1% growth over baseline. Otherwise, check if change is reasonable!

    EXPECT_LE(static_cast<double>(total_size), static_cast<double>(baseline_total_size) * margin)
            << "Total encoder storage size (" << total_size << " bytes) grew by more than 1% over "
            << "the baseline (" << baseline_total_size << " bytes). If this growth is intentional, "
            << "update baseline_total_size to " << total_size << ".";

    EXPECT_LE(static_cast<double>(total_size_819k), static_cast<double>(baseline_819k_total_size) * margin)
            << "Total 819k-encoder storage size (" << total_size_819k << " bytes) grew by more than 1% "
            << "over the baseline (" << baseline_819k_total_size << " bytes). If this growth is "
            << "intentional, update baseline_819k_total_size to " << total_size_819k << ".";

    std::cout << "Total encoder storage: " << total_size << " bytes (baseline: " << baseline_total_size
              << ", +" << (100.0 * static_cast<double>(total_size) / static_cast<double>(baseline_total_size) - 100.0)
              << "%)" << std::endl;
}

TEST(test_prebuilt_codes, basic_example_code_choice_runtime) {
    unsigned seed = 42; // seed for PRNG

    // If the code choice is done at RUNTIME (will usually be the case, e.g. because QBER is known only at runtime),
    // need to provide containers for input and output with correct sizes, otherwise an exception will be thrown.
    // if we don't know the size of `key` at runtime, it will probably be a `std::vector`, not a `std::array`.
    std::size_t code_id = 0;  // Not `constexpr`, let's say we only know this value at runtime!

    std::vector<std::uint8_t> key_vec{};
    key_vec.resize(LDPC4QKD::get_input_size(code_id)); // get size at runtime!
    noise_bitstring_inplace(key_vec, 0.5, seed);  // create a random key

    // allocate a buffer for the syndrome (runtime-known length).
    std::vector<std::uint8_t> syndrome_vec{};
    syndrome_vec.resize(LDPC4QKD::get_output_size(code_id));

    // However,  if `key_vec.size()` and `syndrome_vec.size()` aren't exactly right for the chosen code,
    // then the above will throw `std::out_of_range` exception!
    LDPC4QKD::encode_with(code_id, key_vec, syndrome_vec);

    // Use the non-generic version of `encode_with`:
    std::cout << "Syndrome of runtime known size " << syndrome_vec.size() << std::endl;
    for (auto v: syndrome_vec) {
        std::cout << static_cast<int>(v) << ' ';  // print syndrome bits
    }
    std::cout << std::endl;
}

TEST(test_prebuilt_codes, basic_example_code_choice_runtime_vectorbool) {
    unsigned seed = 42; // seed for PRNG

    // same thing with `vector<bool>`
    std::size_t code_id = 0;  // Not `constexpr`, let's say we only know this value at runtime!

    std::vector<bool> key_vec{};
    key_vec.resize(LDPC4QKD::get_input_size(code_id)); // get size at runtime!
    noise_bitstring_inplace(key_vec, 0.5, seed);  // create a random key

    // allocate a buffer for the syndrome (runtime-known length).
    std::vector<bool> syndrome_vec{};
    syndrome_vec.resize(LDPC4QKD::get_output_size(code_id));

    // However,  if `key_vec.size()` and `syndrome_vec.size()` aren't exactly right for the chosen code,
    // then the above will throw `std::out_of_range` exception!
    LDPC4QKD::encode_with(code_id, key_vec, syndrome_vec);

    // Use the non-generic version of `encode_with`:
    std::cout << "Syndrome of runtime known size " << syndrome_vec.size() << std::endl;
    for (auto v: syndrome_vec) {
        std::cout << static_cast<int>(v) << ' ';  // print syndrome bits
    }
    std::cout << std::endl;
}

TEST(test_prebuilt_codes, basic_example_code_choicecomptime) {
    unsigned seed = 42; // seed for PRNG

    // If the block size and syndrome size are known at compile time, we can use fixed-length buffers (`std::array`)
    // If any of the containers used for key or syndrome has a compile-time known size
    // (e.g. `std::array` or `std::span`), then this method MUST be used!
    std::array<std::uint8_t, AutogenLDPC_QC_2048x4096_0c809c3::N * AutogenLDPC_QC_2048x4096_0c809c3::expansion_factor> key_arr{};
    noise_bitstring_inplace(key_arr, 0.5, seed);  // create a random key

    // allocate a buffer for the syndrome
    std::array<std::uint8_t, AutogenLDPC_QC_2048x4096_0c809c3::M * AutogenLDPC_QC_2048x4096_0c809c3::expansion_factor> syndrome{};

//    encoder1.encode(key_arr, syndrome);  // this would work. It's just using a concrete encoder object.

    // This also works and does the same thing.
    // Use this way when LDPC code choice is known at compile time.
    // That's because `key` and `syndrome` have the correct sizes for `code_id = 0` (which is precisely `encoder1`)
    constexpr int code_id = 3;  // HAS to be `constexpr`!
    LDPC4QKD::encode_with_static<code_id>(key_arr, syndrome);

    // this will not compile because `key` has a compile-time known size.
    // The template instantiation will try to compile against **each** available code (although only one is used).
    // Since the codes have different sizes, this will fail on most codes and give a compiler error.
//    LDPC4QKD::encode_with(code_id, key, syndrome);  // error: no matching function for call to ‘std::span<...>::span(...)

    std::cout << "Syndrome of compile-time known size " << syndrome.size() << std::endl;
    for (auto v: syndrome) {
        std::cout << static_cast<int>(v) << ' ';  // print syndrome bits
    }
    std::cout << std::endl;
}
