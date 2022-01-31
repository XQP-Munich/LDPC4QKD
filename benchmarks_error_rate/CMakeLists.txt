# CMake file for building the error correction demo.

cmake_minimum_required(VERSION 3.19)

project(ErrorCorrectionDemoForQKD)


# ---------------------------------------------------------------------------------- Rate adapted performance simulation
add_executable(rate_adapted_simulation main_rate_adapted_simulation.cpp
        code_simulation_helpers.hpp)

target_compile_features(rate_adapted_simulation PUBLIC cxx_std_17)

target_link_libraries(rate_adapted_simulation
        PRIVATE
        # build options
        compiler_warnings
        project_options

        # libraries
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
        # build options
        compiler_warnings
        project_options

        # libraries
        LDPC4QKD::LDPC4QKD
        )

target_include_directories(rate_adapted_fer
        PRIVATE
        )
