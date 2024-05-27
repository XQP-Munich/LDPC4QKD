//
// Created by Adomas Baliuka on 02.05.24.
//

// Google Test framework
#include <gtest/gtest.h>
#include "helpers_for_testing.hpp"

// Standard library
#include <iostream>

// To be tested
#include "encoder_advanced.hpp"

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

TEST(test_encoder_advanced, basic_example_code_choice_runtime) {
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

TEST(test_encoder_advanced, basic_example_code_choice_runtime_vectorbool) {
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

TEST(test_encoder_advanced, basic_example_code_choicecomptime) {
    unsigned seed = 42; // seed for PRNG

    // If the block size and syndrome size are known at compile time, we can use fixed-length buffers (`std::array`)
    // If any of the containers used for key or syndrome has a compile-time known size
    // (e.g. `std::array` or `std::span`), then this method MUST be used!
    std::array<std::uint8_t, AutogenLDPC_QC::N * AutogenLDPC_QC::expansion_factor> key_arr{};
    noise_bitstring_inplace(key_arr, 0.5, seed);  // create a random key

    // allocate a buffer for the syndrome
    std::array<std::uint8_t, AutogenLDPC_QC::M * AutogenLDPC_QC::expansion_factor> syndrome{};

//    encoder1.encode(key_arr, syndrome);  // this would work. It's just using a concrete encoder object.

    // This also works and does the same thing.
    // Use this way when LDPC code choice is known at compile time.
    // That's because `key` and `syndrome` have the correct sizes for `code_id = 0` (which is precisely `encoder1`)
    constexpr int code_id = 0;  // HAS to be `constexpr`!
    LDPC4QKD::encode_with<code_id>(key_arr, syndrome);

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
