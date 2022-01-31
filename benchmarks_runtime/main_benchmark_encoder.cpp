//
// Created by alice on 08.06.21.
//

// Standard library
#include <iostream>
#include <random>

// Google Benchmark library
#include <benchmark/benchmark.h>

// Project scope
#include "encoder.hpp"


template<typename T, std::size_t N>
void noise_bitstring_inplace(std::array<T, N> &src, double err_prob, unsigned int seed = 0) {
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


std::array<bool, AutogenLDPC::N> input{};
std::array<bool, AutogenLDPC::M> output{};
std::array<bool, AutogenLDPC::M/2 + 100> ra_output{};

static void BM_encode_benchmark_no_rate_adaption(benchmark::State& state) {
    // Perform setup here

    for (auto _ : state) {
        state.PauseTiming();
        noise_bitstring_inplace(input, 0.5, static_cast<unsigned int>(state.range(0)));
        state.ResumeTiming();

        benchmark::DoNotOptimize(input);
        benchmark::DoNotOptimize(output);

        // This code gets timed
        LDPC4QKD::encode(input, output);

        benchmark::ClobberMemory();
    }
}
BENCHMARK(BM_encode_benchmark_no_rate_adaption)->Unit(benchmark::kMicrosecond)->Range(8, 8<<10);


static void BM_encode_benchmark_with_rate_adaption(benchmark::State& state) {
    // Perform setup here

    for (auto _ : state) {
        state.PauseTiming();
        noise_bitstring_inplace(input, 0.5, static_cast<unsigned int>(state.range(0)));
        state.ResumeTiming();

        benchmark::DoNotOptimize(input);
        benchmark::DoNotOptimize(output);
        benchmark::DoNotOptimize(ra_output);

        LDPC4QKD::encode(input, output);
        LDPC4QKD::rate_adapt(output, ra_output);

        benchmark::ClobberMemory();
    }
}
BENCHMARK(BM_encode_benchmark_with_rate_adaption)->Unit(benchmark::kMicrosecond)->Range(8, 8<<10);


static void BM_only_rate_adaption(benchmark::State& state) {
    // Perform setup here

    for (auto _ : state) {
        state.PauseTiming();
        noise_bitstring_inplace(output, 0.5, static_cast<unsigned int>(state.range(0)));
        state.ResumeTiming();

        benchmark::DoNotOptimize(output);
        benchmark::DoNotOptimize(ra_output);

        LDPC4QKD::rate_adapt(output, ra_output);

        benchmark::ClobberMemory();
    }
}
BENCHMARK(BM_only_rate_adaption)->Unit(benchmark::kMicrosecond)->Range(8, 8<<10);


// Run the benchmark
BENCHMARK_MAIN();
