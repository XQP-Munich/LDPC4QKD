[![workflow](https://github.com/XQP-Munich/LDPC4QKD/actions/workflows/ci-cmake_tests.yml/badge.svg)](https://github.com/XQP-Munich/LDPC4QKD/actions)
[![codecov](https://codecov.io/gh/XQP-Munich/LDPC4QKD/branch/main/graph/badge.svg?token=GV9453ZM42)](https://codecov.io/gh/XQP-Munich/LDPC4QKD)
[![License](https://img.shields.io/github/license/XQP-Munich/LDPC4QKD)](./LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5579246.svg)](https://doi.org/10.5281/zenodo.5579246)
# LDPC4QKD: LDPC codes for rate adaptive distributed source coding 

## Overview

This repository aims to solve the following problem: Suppose Alice and Bob each have a bit-string of length _N_.
Bobs bit-string is the result of randomly flipping each bit of Alice with probability _p_.
I.e., it is the output of a [binary symmetric channel](https://en.wikipedia.org/wiki/Binary_symmetric_channel) with **channel parameter** _p_. 
(Note: similar considerations apply for soft inputs instead of bit-strings and a [Gaussian channel](https://en.wikipedia.org/wiki/Additive_white_Gaussian_noise).)
The goal is for Alice to send (via a noise-less channel) to Bob a message (called syndrome), such that Bob can recover Alice's bit-string given his bit-string and the syndrome received from Alice.

There are two important metrics, which we want to minimize: 

- the probability that Bob will fail to recover Alice's bit-string, which is called the frame error rate (FER)
- the length of the syndrome

We define the syndrome to be the matrix-vector product (modulo 2) of a sparse binary matrix (LDPC matrix) and Alice's bit-string. 
In this case, Bob can use an inference algorithm (loopy belief propagation) to obtain a guess of Alice's bit-string. If the LDPC matrix has size _M_ x _N_, we call _M_ / _N_ the **rate** of the LDPC matrix and _N_ the **block size**. 
(Note: the term rate is used differently in forward error correction, where it means 1 - _M_ / _N_.) Minimizing the length of the syndrome actually means minimizing the rate for a given channel parameter.
Note that, for any given channel parameter, there is a trade-off between small rate, small FER and small block size.
Furthermore, there are theoretical limits as to how small the rate can be (the Slepian Wolf limit, sometimes called Shannon limit, applies to assymptotically large bloc size.
Limits for finite block sizes also exist but are more complicated).

### Contrast with forward error correction

Contrary to how forward error correction works, distributed source coding does not use generator matrices for the LDPC code. 
This repository does not provide generator matrices (though calculating them from the parity check matrices is straightforward).

Our decoder implementation operates on a bit-string (noisy version of true message) and its correct syndrome. 
This is different from what is used for forward error correction (as in e.g. [AFF3CT](https://github.com/aff3ct/aff3ct)), where the decoder operates on only the noisy codeword.
(The noisy codeword is the result of transmitting the codeword (true message encoded using a generator matrix) via a noisy channel.)


## How to use
In order to use all the functionality provided in this repository, install
- CMake (at least version 3.19)
- C++ compiler (supporting C++20; parts of the project also work with only C++17)
- Julia (at least version 1.6)

For how to use the provided Julia code, see directory `codes`. No familiarity with the Julia programming language is required.

To build the C++ executables (except for unit tests) using CMake (Julia not required), execute

        git clone <github url>
        cd LDPC4QKD
        cmake -S . -B build -DBUILD_UNIT_TESTS=OFF
        cmake --build build --config Release
        
All executables will be built inside the `build` folder.

For a demo of how to use the C++ header-only library, see the `examples` directory, which shows a basic example ("demo_error_correction") of how to use the C++ header containing the decoder.
This example is built by CMake (executable `build/examples/demo_error_correction`.
Note: the executable produces no output; look at the C++ source code to see how the header can be used.

For more advanced examples, looking at the unit tests may be helpful.

## How to contribute
This repository is actively maintained. 
Issues and pull requests will be responded to and processed.

Let us know if you're having problems with the provided materials, wish to contribute new ideas or need modifications of our work for your application.

## List of contents

### Data files
- A number of LDPC codes (a list and the actual LDPC matrices are in the folder `codes`).
  Their parity check matrices are stored in a custom file format (called `qccsc.json`).

- Simulations results done using [AFF3CT](https://github.com/aff3ct/aff3ct), showing FER of the LDPC matrices at various channel parameters.

- Simulation results done using the decoder in this repository, showing FER of LDPC matrices, their rate adapted versions, and average rate under rate adaption (Work in progress!).

- For each LDPC matrix, a specification of rate adaption. 
  This is a list of pairs of row indices of the matrix that are combined (added mod 2) in each rate adaption step.

### Julia code
- Contained in the folder `codes`.

- Enables loading the LDPC matrices from files.
  This uses our custom Julia library [LDPCStorage.jl](https://github.com/XQP-Munich/LDPCStorage.jl)) (installed automatically by the Julia package manager).

- Enables saving the compressed sparse column (CSC) representation (both binary and quasi-cyclic).
  Also supports exporting to other standard formats, such as `alist`.

### C++ code
- Basic LDPC decoder using belief propagation (BP).
  Contained in a single header file (`src/rate_adaptive_code.hpp`) and easy to use.
  Can perform syndrome computation as well. 
  The LDPC code can be embedded into the executable or loaded from a file at program runtime.

- Utility functions for reading LDPC matrices (from `.cscmat` files) and rate adaption (from `.csv` files).
  These functions are contained in `src/read_cscmat_format.hpp`.
  Note that these functions require the fully explicit binary form of the LDPC matrix.
  The LDPC codes given inside `codes/ldpc` use an even more memory-efficient storage.

- LDPC matrices can be stored within the executable.
  There are two ways:
  - include as a header file defining constant data (for example `tests/fortest_autogen_ldpc_matrix_csc.hpp`).
    Julia code (see folder `codes`) is provided to generate such a C++ header file for any LDPC matrix.
  - C++ headers storing QC-LDPC matrices in terms of their quasi-cyclic exponents.
    This is very efficient in terms of binary size.
    See `src/autogen_ldpc_QC.hpp` (partially auto-generated using the Julia code).
    Note: this part is new and may still change significantly in future versions.
    It also requires C++20 and `src/encoder_advanced.hpp`.

- For applications that only require syndrome computation but no decoding, we provide a separate implementation for multiplication of a sparse binary matrix and a dense binary vector (LDPC syndrome computation).
  See `src/encoder.hpp` (old) or `src/encoder_advanced.hpp` (new).
  **This is a very specific application, which you probably don't care about initially.**

## Planned features and improvements

If you need some feature for your applications, let us know, e.g. by creating an issue.

- [ ] Reports and benchmarks on LDPC code quality and decoding performance
  + [x] Simulation programs with command line interfaces
  + [x] Frame error rate simulations for rate adapted codes (including special case of no rate adaption)
  + [x] Critical rate (codeword-averaged minimum leak rate for successful decoding) computation for rate adapted codes
  + [ ] Automatic performance reports with code that generates plots
- [ ] LDPC codes
  + [x] 3 LDPC codes each (different block sizes) for leak rates 1/2 and 1/3
  + [ ] More sizes, more rates (coming soon!)
- [ ] Decoding and decoding algorithms
  + [x] Basic belief propagation (BP) decoder for Slepian-Wolf setting
  + [ ] Decoder performance improvements (look at [AFF3CT](https://github.com/aff3ct/aff3ct) for inspiration), plausibly achieve 2x runtime speedup at same decoding accuracy
  + [ ] Decoding on GPU
  + [x] Encoder that can be used separately from encoder (e.g. for embedded applications)
  + [x] Save memory by storing QC-exponents of structured codes, rather than CSC storage
  + [x] Encoder that directly uses CSC-storage of QC-exponents instead of expanding to CSC storage of binary matrix

## Attributions

Some of the simulation/benchmarking programs use
- [CmdParser](https://github.com/FlorianRappl/CmdParser), a simple command line argument parser (MIT license, the sources are included in this repository at `external/CmdParser-91aaa61e`). Copyright (c) 2015 - 2016 Florian Rappl
- [json](https://github.com/nlohmann/json), a JSON parser (MIT License and others, see header `external/json-6af826d/json.hpp`). © 2013-2022 Niels Lohmann.
