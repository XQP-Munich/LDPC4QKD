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

# building unit tests. Uses Google Tests framework, which has to be installed on the target system.
# NOTE: For the file-IO tests to work, `tests` must be the cwd of the test execution. If this fails, either adapt the
# path variable tests_folder_path in test_main.cpp or change the cwd of the executable.
add_executable(unit_tests_error_correction test_main.cpp
        # -------- Actual Unit tests --------
        test_encoder.cpp
        test_rate_adaptive_code.cpp
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
