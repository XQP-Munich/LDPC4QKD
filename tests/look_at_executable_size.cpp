//
// Created by alice on 26.04.21.
// This should be built as a stand-alone application.
// It's only purpose is to see how much memory it takes to save the LDPC code within the executable.
//

#include <iostream>
#include <array>
#include <cstdint>
#include <random>

#include "fortest_autogen_ldpc_matrix_csc.hpp"

using namespace AutogenLDPC;

template<typename Bit>
constexpr void encode(const std::array<Bit, N>& in, std::array<Bit, M>& out) {

    for (std::size_t col = 0; col < N; col++) {
        for (std::size_t j = colptr[col]; j < colptr[col + 1]; j++) {
            out[row_idx[j]] = (static_cast<bool>(out[row_idx[j]]) != static_cast<bool>(in[col]));
        }
    }
}

template<typename T, std::size_t N>
void noise_bitstring_inplace(std::array<T, N> &src, double err_prob, unsigned int seed = 0) {
    std::mt19937_64 rng{seed}; // hard-coded seed for testing purposes.

    std::bernoulli_distribution distribution(err_prob);

    for (std::size_t i = 0; i < N; i++) {
        if (distribution(rng)) {
            src[i] = !src[i];
        } else {
            src[i] = src[i];
        }
    }
}

template<class T>
void print_nz_inds(const T &vec) {
    std::size_t i{};
    for (auto val : vec) {
        if (val != 0) {
            std::cout << i << ' ';
        }
        i++;
    }
    std::cout << std::endl;
}

using Bit = bool;
std::array<Bit, N> input{};
std::array<Bit, M> output{};

int main() {
    noise_bitstring_inplace(input, 0.5);
    encode(input, output);

    print_nz_inds(output);
}
