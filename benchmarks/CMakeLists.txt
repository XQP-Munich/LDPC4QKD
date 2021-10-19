# CMake file for building the error correction demo.

cmake_minimum_required(VERSION 3.19)

project(ErrorCorrectionDemoForQKD)

# google benchmark library. Only required for the targets called "benchmark", i.e., for runtime speed benchmarking.
find_package(benchmark REQUIRED)

# --------------------------------------------------------------------------------------------------- Encoder Benchmarks
add_executable(benchmark_encoder main_benchmark_encoder.cpp
        )

target_compile_features(benchmark_encoder PUBLIC cxx_std_17)

target_link_libraries(benchmark_encoder
        PRIVATE
        LDPC4QKD::LDPC4QKD
        benchmark::benchmark
        )

target_include_directories(benchmark_encoder
        PRIVATE
        )

# --------------------------------------------------------------------------------------------------- Decoder Benchmarks
add_executable(benchmark_decoder main_benchmark_decoder.cpp
        )

target_compile_features(benchmark_decoder PUBLIC cxx_std_17)

target_link_libraries(benchmark_decoder
        PRIVATE
        LDPC4QKD::LDPC4QKD
        benchmark::benchmark
        )

target_include_directories(benchmark_decoder
        PRIVATE
        )

# --------------------------------------------------------------------------------------------- Rate Adaption Benchmarks
add_executable(benchmark_ra main_benchmark_ra.cpp
        )

target_compile_features(benchmark_ra PUBLIC cxx_std_17)

target_link_libraries(benchmark_ra
        PRIVATE
        LDPC4QKD::LDPC4QKD
        benchmark::benchmark
        )

target_include_directories(benchmark_ra
        PRIVATE
        )

# ---------------------------------------------------------------------------------- Rate adapted performance simulation
add_executable(rate_adapted_simulation main_rate_adapted_simulation.cpp
        code_simulation_helpers.hpp)

target_compile_features(rate_adapted_simulation PUBLIC cxx_std_17)

target_link_libraries(rate_adapted_simulation
        PRIVATE
        LDPC4QKD::LDPC4QKD
        )

target_include_directories(rate_adapted_simulation
        PRIVATE
        )

# ------------------------------------------------------------- frame error rate (FER) simulation (allows rate adaption)
add_executable(rate_adapted_fer main_rate_adapted_fer.cpp
        code_simulation_helpers.hpp)

target_compile_features(rate_adapted_fer PUBLIC cxx_std_17)

target_link_libraries(rate_adapted_fer
        PRIVATE
        LDPC4QKD::LDPC4QKD
        )

target_include_directories(rate_adapted_fer
        PRIVATE
        )
