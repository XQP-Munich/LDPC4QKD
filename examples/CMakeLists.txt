# CMake file for building the error correction demo.

cmake_minimum_required(VERSION 3.19)

project(ErrorCorrectionDemo)

# ------------------------------------------------- BASIC DEMO ---------------------------------------------------------
add_executable(demo_error_correction main_demo_error_correction.cpp
        )

target_compile_features(demo_error_correction PUBLIC cxx_std_17)

target_link_libraries(demo_error_correction
        PRIVATE
        LDPC4QKD::LDPC4QKD
        )
