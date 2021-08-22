[![workflow](https://github.com/XQP-Munich/LDPC4QKD/actions/workflows/ci-cmake_tests.yml/badge.svg)](https://github.com/XQP-Munich/LDPC4QKD/actions)
[![License](https://img.shields.io/github/license/XQP-Munich/LDPC4QKD)](./LICENSE)

# LDPC4QKD: Rate adaptive LDPC-based methods for distributed source coding.

Note: This repository is still missing some features which will be added very soon.

## Overview

The problem solved by distributed source coding is the following. Suppose Alice and Bob 

Contrary to how forward error correction works, generator matrices for the LDPC code are not used at all for distributed source coding. This repository does not provide generator matrices (although calculating them from the parity check matrix is straightforward, s)
Our decoder operates on a bit-string (noisy version of true message) and its correct syndrome. This is slightly different from what is used for forward error correction (as in e.g. AFF3CT), where the decoder operates on only the noisy codeword (which is the result transmitting the codeword, which in turn is the encoding the true message using a generator matrix).

## How to contribute
This repository is actively maintained. Issues and pull requests will be responded to and processed.

If you're having problems with the provided materials, want to contribute new ideas or need modifications or our work for your application, please contact us.

## List of contents

### Data files
- A number of LDPC codes. Their parity check matrices are stored in a custom file format (called `CSCMAT`).
- Simulations results (done using AFF3CT) showing the quality of the LDPC matrices.
- For each LDPC matrix, a specification of rate adaption. This is a list of pairs of row indices of the matrix that are combined (added mod 2) in each rate adaption step. 

### Julia code to
- load the LDPC matrices from `CSCMAT` files
- save the compressed sparse column (CSC) representation and also to export to other standard formats, such as `alist`.

### C++ code
- Basic LDPC decoder using belief propagation (BP). Contained in a single header file (`PATH`) and easy to use. Can perform syndrome computation as well. The LDPC code can be embedded into the executable or loaded from a file at program runtime.
- For applications that only require syndrome computation but no decoding, we provide a separate implementation for multiplication of a sparse binary matrix and a dense binary vector (LDPC syndrome computation). This is also contained in a single header file (`PATH`).
- For this purpose, the LDPC matrix can be stored within the executable. Julia code is provided to automatically generate header files storing an LDPC code in constant arrays.

TODO add paths

## How to use

Start from the `examples` directory, which shows how to use the C++ header containing the decoder. 
TODO


## Planned features and improvements
- Code to automatically generate reports on quality of all LDPC codes and rate adapted performance (work in progress)
- More and better LDPC codes with more sizes and rates
- The BP implementation is not state-of-the-art or heavily optimized. We may improve it in the future.
- Using compressed sparse row format may have advantages over the current compressed sparse column format. This needs to be tested.
- the quasi-cyclic structure of LDPC matrices can be exploited to store the matrices more efficiently. For the encoder, this may be combinable with storing individual bits as integers rather than booleans, again saving some memory. The encoding in this case would consist of adding bit-shifted integers.


## Appendix

### Table of LDPC matrices

TODO