cmake_minimum_required(VERSION 3.19)

project(LDPC4QKDLibrary)

add_library(LDPC4QKD INTERFACE
        LDPC4QKD/rate_adaptive_code.hpp  # contains decoder and encoder. Does not need any other files to work.
        LDPC4QKD/autogen_ldpc_QC.hpp
        LDPC4QKD/encoder.hpp # contains only minimal encoder that needs static storage of the LDPC matrix
        LDPC4QKD/encoder_advanced.hpp # REQUIRES C++20!!! advanced encoder (QC-enabled and constexpr objects). Handles storage.
        LDPC4QKD/read_ldpc_file_formats.hpp # helper methods to generate static storage (not needed to use encoder/decoder class).
)

target_compile_features(LDPC4QKD INTERFACE cxx_std_20)

target_include_directories(LDPC4QKD
        INTERFACE
        ${CMAKE_CURRENT_LIST_DIR}
        )

# Create alias with "namespace". This is the target to be used to link to the library in higher level projects.
add_library(LDPC4QKD::LDPC4QKD ALIAS LDPC4QKD)
