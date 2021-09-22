//
// Created by alice on 20.09.21.
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
#include "autogen_rate_adaption.hpp"


LDPC4QKD::RateAdaptiveCode<bool> get_code_big_wra() {
    std::vector<std::uint32_t> colptr(AutogenLDPC::colptr.begin(), AutogenLDPC::colptr.end());
    std::vector<std::uint16_t> row_idx(AutogenLDPC::row_idx.begin(), AutogenLDPC::row_idx.end());
    std::vector<std::uint16_t> rows_to_combine(AutogenRateAdapt::rows.begin(), AutogenRateAdapt::rows.end());
    return LDPC4QKD::RateAdaptiveCode<bool>(colptr, row_idx, rows_to_combine);
}


static void BM_decode_benchmark_set_rate(benchmark::State& state) {
    for (auto _ : state) {
        state.PauseTiming();
        auto H = get_code_big_wra();
        state.ResumeTiming();

        benchmark::DoNotOptimize(H.getPosVarn());
        H.set_rate(state.range(0));
        benchmark::ClobberMemory();
    }
}
BENCHMARK(BM_decode_benchmark_set_rate)->Unit(benchmark::kMicrosecond)->Range(0, 1024);


// Run the benchmark
BENCHMARK_MAIN();
