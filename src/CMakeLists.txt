cmake_minimum_required(VERSION 3.16)

project(LDPC4QKDLibrary)


if(NOT EXISTS "${CMAKE_CURRENT_LIST_DIR}/autogen_ldpc_matrix_csc.hpp")
    message(FATAL_ERROR "C++ code containing the LDPC matrix not found.
    Use the Julia scripts to automatically generate it form an .alist file!")
endif(NOT EXISTS "${CMAKE_CURRENT_LIST_DIR}/autogen_ldpc_matrix_csc.hpp")


if(NOT EXISTS "${CMAKE_CURRENT_LIST_DIR}/autogen_rate_adaption.hpp")
    message(FATAL_ERROR "C++ code containing the rate adaption for the LDPC matrix not found.
    Use Julia to automatically generate it form an appropriate CSV file!")
endif(NOT EXISTS "${CMAKE_CURRENT_LIST_DIR}/autogen_rate_adaption.hpp")


add_library(LDPC4QKD INTERFACE
        rate_adaptive_code.hpp
        encoder.hpp
        )

target_compile_features(LDPC4QKD INTERFACE cxx_std_17)

target_include_directories(LDPC4QKD
        INTERFACE
        ${CMAKE_CURRENT_LIST_DIR}
        )

# Create alias with "namespace". This is the target to be used to link to the library in higher level projects.
add_library(LDPC4QKD::LDPC4QKD ALIAS LDPC4QKD)
