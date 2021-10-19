#
# CMakeLists.txt file for building unit tests.
# Adapted from https://github.com/pabloariasal/modern-cmake-sample/blob/master/libjsonutils/test/CMakeLists.txt
#

cmake_minimum_required(VERSION 3.19)


include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG        release-1.11.0
)
FetchContent_MakeAvailable(googletest)


# building unit tests. Uses Google Tests framework, which has to be installed on the target system.
# NOTE: For the file-IO tests to work, `tests` must be the cwd of the test execution. If this fails, either adapt the
# path variable tests_folder_path in test_main.cpp or change the cwd of the executable.
add_executable(unit_tests_error_correction test_main.cpp
        # -------- Actual Unit tests --------
        test_encoder.cpp
        test_rate_adaptive_code.cpp
        test_read_scmat_format.cpp

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


get_filename_component(LDPC4QKD_TESTS_DIR ${CMAKE_CURRENT_LIST_DIR} DIRECTORY)
get_filename_component(LDPC4QKD_SRC_DIR "${LDPC4QKD_TESTS_DIR}/src" ABSOLUTE)
# Another executable that is only used to see what size the compiled binary has for a given LDPC code,
# when embedding the code into the binary.
add_executable(look_at_executable_size look_at_executable_size.cpp)
target_include_directories(look_at_executable_size
        PRIVATE
        "${LDPC4QKD_SRC_DIR}/"
        )

# Copy the .cscmat file into the directory containing the tests binary.
# This is required to test the .cscmat reader code.
get_filename_component(CSCMAT_TEST_FILE_PATH
        ${CMAKE_CURRENT_SOURCE_DIR}/LDPC_code_for_testing_2048x6144.cscmat
        REALPATH
        )
message(STATUS "Copying file ${CSCMAT_TEST_FILE_PATH} into directory ${CMAKE_CURRENT_BINARY_DIR}.")
file(COPY ${CSCMAT_TEST_FILE_PATH}
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
        )
