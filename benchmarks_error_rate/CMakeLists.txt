# CMake file for building the error correction demo.

cmake_minimum_required(VERSION 3.19)

project(ErrorCorrectionDemoForQKD)


# ---------------------------------------------------------------------------------- Rate adapted performance simulation
add_executable(critical_rate_simulation main_critical_rate_simulation.cpp
        code_simulation_helpers.hpp)

target_compile_features(critical_rate_simulation PUBLIC cxx_std_17)

target_link_libraries(critical_rate_simulation
        PRIVATE
        # build options
        compiler_warnings
        project_options

        # libraries
        LDPC4QKD::LDPC4QKD
        )

target_include_directories(critical_rate_simulation
        PRIVATE
        "${CMAKE_CURRENT_LIST_DIR}/../"
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
        "${CMAKE_CURRENT_LIST_DIR}/../"
        )
