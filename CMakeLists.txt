# CMake file building demos and tests for the error correction library.

cmake_minimum_required(VERSION 3.19)
message("CMake version: ${CMAKE_VERSION}")
project(LDPC4QKD
        VERSION 0.0.0 # This version is meaningless at the moment! TODO generate Cmake-Project-version dynamically
        LANGUAGES C CXX
        DESCRIPTION "Rate adaptive distributed source coding using LDPC codes."
        )

# Allows `include`'ing custom CMake modules from this directory.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")


# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------- CMAKE OPTIONS --------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

# !!!!!
# Note:
# CMake saves options in the cache (`CMakeCache.txt`).
# If CMake `option` default values are changed in this `CMakeLists.txt` file, the old values from the cache will
# continue to be used. Delete the cache for such changes to take effect.
# !!!!!

option(DEBUG_MESSAGES_ENABLED
        "Enables the output debug messages (compile definitions using macros in the C++ code)
        "
        ON
        )

option(BUILD_RUNTIME_BENCHMARKS
        "Enables building of runtime performance benchmarks.
        The benchmarks use the Google Benchmark framework, which must be installed.
        "
        OFF
        )

option(BUILD_UNIT_TESTS
        "Enables building unit tests.
        The unit tests use the Google Test framework, which is downloaded and installed automatically!
        "
        ON
        )

if (DEBUG_MESSAGES_ENABLED)
    message(STATUS "#define compile definition DEBUG_MESSAGES_ENABLED.")
    add_compile_definitions(DEBUG_MESSAGES_ENABLED)
endif ()



# Interface library (without sources) that carries compiler/build flags.
add_library(project_options INTERFACE)
add_library(compiler_warnings INTERFACE)

include(BuildConfig)  # custom CMake module (see folder `cmake/Modules`)
enable_sanitizers(project_options) # adds coverage if specified by flag ENABLE_COVERAGE

include(CompilerWarnings)
set_project_warnings(compiler_warnings)

# ----------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------- MISCELLANEOUS ---------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

# This will allow us to print a feature summary.
# https://cmake.org/cmake/help/v3.11/module/FeatureSummary.html
include(FeatureSummary)


# Required to enable testing for link-time optimization.
# https://cmake.org/cmake/help/v3.11/module/CheckIPOSupported.html
include(CheckIPOSupported)


# Use link-time optimization if allowed.
check_ipo_supported(RESULT ipo_supported)
if (ipo_supported)
    set_property(GLOBAL PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    message(STATUS "Link-time optimization (IPO) is enabled.")
else ()
    message(WARNING "Link-time optimization is not supported.")
endif ()


# We attempt to use ccache to speed up the build.
find_program(CCACHE_FOUND "ccache")
if (CCACHE_FOUND)
    message(STATUS "Using ccache for building.")
    set_property(GLOBAL PROPERTY LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
    message(STATUS "ccache is being used to speed up the build.")
else ()
    message(WARNING "Using ccache to speed up the build is not possible.")
endif (CCACHE_FOUND)

# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------- INCLUDE SUBDIRECTORIES -----------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------

add_subdirectory(src)

add_subdirectory(examples)

if (BUILD_UNIT_TESTS)
    # Enable testing. This is not necessary for using GTest,
    # only provides integration with CTest (CMake test suite)
    # Note that this line has to be located before adding the sub-directory with the tests.
    # Having this line in the test/CMakeLists.txt file doesn't work (I don't know why).
    enable_testing()

    add_subdirectory(tests)
else (BUILD_UNIT_TESTS)
    message(STATUS "Unit tests are not built, as specified by user.")
endif (BUILD_UNIT_TESTS)


if (BUILD_RUNTIME_BENCHMARKS)
    add_subdirectory(benchmarks)
else (BUILD_RUNTIME_BENCHMARKS)
    message(STATUS "Performance benchmarks are not built, as specified by user.")
endif (BUILD_RUNTIME_BENCHMARKS)
