# CMake file for building the error correction demo.

cmake_minimum_required(VERSION 3.19)

project(ErrorCorrectionDemoForQKD)

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

# ------------------------------------------------------------------------------------------ Frame error rate simulation
add_executable(frame_error_rate_simulation main_frame_error_rate_simulation.cpp
        code_simulation_helpers.hpp)

target_compile_features(frame_error_rate_simulation PUBLIC cxx_std_17)

target_link_libraries(frame_error_rate_simulation
        PRIVATE
        LDPC4QKD::LDPC4QKD
        )

target_include_directories(frame_error_rate_simulation
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
