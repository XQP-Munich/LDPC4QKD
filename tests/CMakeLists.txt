#
# CMakeLists.txt file for building unit tests.
# Adapted from https://github.com/pabloariasal/modern-cmake-sample/blob/master/libjsonutils/test/CMakeLists.txt
#

cmake_minimum_required(VERSION 3.17)


include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        release-1.11.0
)
FetchContent_MakeAvailable(googletest)


if(NOT EXISTS "${CMAKE_CURRENT_LIST_DIR}/fortest_autogen_ldpc_matrix_csc.hpp")
    message(WARNING "C++ code containing the LDPC matrix for unit tests not found.
    Use the Julia scripts to automatically generate it form an .alist file! A minimal example will be generated automatically.")

    add_custom_command(
            OUTPUT fortest_autogen_ldpc_matrix_csc.hpp
            DEPENDS "${CMAKE_SOURCE_DIR}/src/autogen_ldpc_matrix_csc_example.hpp"
            COMMAND ${CMAKE_COMMAND} -E copy
            "${CMAKE_SOURCE_DIR}/src/autogen_ldpc_matrix_csc_example.hpp"
            "${CMAKE_CURRENT_LIST_DIR}/fortest_autogen_ldpc_matrix_csc.hpp"
            VERBATIM
    )
endif(NOT EXISTS "${CMAKE_CURRENT_LIST_DIR}/fortest_autogen_ldpc_matrix_csc.hpp")


if(NOT EXISTS "${CMAKE_CURRENT_LIST_DIR}/fortest_autogen_rate_adaption.hpp")
    message(WARNING "C++ code containing the rate adaption for unit tests not found.
    Use Julia to automatically generate it form an appropriate CSV file! A minimal example will be generated automatically.")

    add_custom_command(
            OUTPUT fortest_autogen_rate_adaption.hpp
            DEPENDS "${CMAKE_SOURCE_DIR}/src/autogen_rate_adaption_example.hpp"
            COMMAND ${CMAKE_COMMAND} -E copy
            "${CMAKE_SOURCE_DIR}/src/autogen_rate_adaption_example.hpp"
            "${CMAKE_CURRENT_LIST_DIR}/fortest_autogen_rate_adaption.hpp"
            VERBATIM
    )
endif(NOT EXISTS "${CMAKE_CURRENT_LIST_DIR}/fortest_autogen_rate_adaption.hpp")


# building unit tests. Uses Google Tests framework, which has to be installed on the target system.
# NOTE: For the file-IO tests to work, `tests` must be the cwd of the test execution. If this fails, either adapt the
# path variable tests_folder_path in test_main.cpp or change the cwd of the executable.
add_executable(unit_tests_error_correction test_main.cpp
        # -------- Actual Unit tests --------
        test_encoder.cpp
        test_rate_adaptive_code.cpp

        # Static data LDPC code used for tests:
        fortest_autogen_ldpc_matrix_csc.hpp
        fortest_autogen_rate_adaption.hpp
        )

target_compile_features(unit_tests_error_correction
        PUBLIC cxx_std_17
        )

target_include_directories(unit_tests_error_correction
        PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}"
        )

target_link_libraries(unit_tests_error_correction
        PRIVATE
        gtest_main
        pthread

        # To be tested:
        LDPC4QKD::LDPC4QKD
        )


# adds tests via CTest
include(GoogleTest)
gtest_discover_tests(unit_tests_error_correction)


# Another executable that is only used to see what size the compiled binary has for a given LDPC code,
# when embedding the code into the binary.
add_executable(look_at_executable_size look_at_executable_size.cpp)
target_include_directories(look_at_executable_size
        PRIVATE
        "${CMAKE_SOURCE_DIR}/src/"
        )
