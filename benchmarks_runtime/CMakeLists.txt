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
        # build options
        compiler_warnings
        project_options

        # libraries
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
        # build options
        compiler_warnings
        project_options

        # libraries
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
        # build options
        compiler_warnings
        project_options

        # libraries
        LDPC4QKD::LDPC4QKD
        benchmark::benchmark
        )

target_include_directories(benchmark_ra
        PRIVATE
        )

if (BUILD_UNIT_TESTS)
    add_test(test_benchmark_ra benchmark_ra)
    add_test(test_benchmark_decoder benchmark_decoder)
    add_test(test_benchmark_encoder benchmark_encoder)
endif (BUILD_UNIT_TESTS)
