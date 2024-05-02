//
// Created by alice on 08.06.21.
//

// Standard library
#include <iostream>
#include <random>

// Google Benchmark library
#include <benchmark/benchmark.h>

// Project scope
#include "rate_adaptive_code.hpp"

// Test cases test against constants known to be correct for the LDPC-matrix defined here:
#include "autogen_ldpc_matrix_csc.hpp"


template<typename T>
void noise_bitstring_inplace(std::vector<T> &src, double err_prob, unsigned int seed = 0) {
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


LDPC4QKD::RateAdaptiveCode get_code_big_nora() {
    std::vector<std::uint32_t> colptr(AutogenLDPC::colptr.begin(), AutogenLDPC::colptr.end());
    std::vector<std::uint16_t> row_idx(AutogenLDPC::row_idx.begin(), AutogenLDPC::row_idx.end());
    return LDPC4QKD::RateAdaptiveCode(colptr, row_idx);
}


static void BM_decode_benchmark_no_rate_adaption(benchmark::State& state) {
    constexpr double bsc_err_p = 0.03; // channel error probability

    // Perform setup here
    auto H = get_code_big_nora();

    std::vector<bool> true_codeword(H.getNCols());
    std::vector<bool> predicted_codeword(H.getNCols());
    std::vector<bool> syndrome(H.get_n_rows_mother_matrix());

    for (auto _ : state) {
        state.PauseTiming();

        noise_bitstring_inplace(true_codeword, 0.5, static_cast<unsigned int>(state.range(0)));
        H.encode_no_ra(true_codeword, syndrome);
        std::vector<bool> x_noised = true_codeword; // distorted data
        noise_bitstring_inplace(x_noised, bsc_err_p);

        const double vlog = log((1 - bsc_err_p) / bsc_err_p);
        std::vector<double> llrs(true_codeword.size());
        for (std::size_t i{}; i < llrs.size(); ++i) {
            llrs[i] = vlog * (1 - 2 * x_noised[i]); // log likelihood ratios
        }
        bool success = false;

        state.ResumeTiming();

        benchmark::DoNotOptimize(predicted_codeword);
        benchmark::DoNotOptimize(syndrome);
        benchmark::DoNotOptimize(llrs);
        benchmark::DoNotOptimize(success);

        // This code gets timed
        success = H.decode_at_current_rate(llrs, syndrome, predicted_codeword);
        if (!success) {
            std::cout << "DECODER FAILED!!!!";
        }

        benchmark::ClobberMemory();
    }
}
BENCHMARK(BM_decode_benchmark_no_rate_adaption)->Unit(benchmark::kMicrosecond)->Range(0, 1<<16);


// Run the benchmark
BENCHMARK_MAIN();
