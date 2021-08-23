[![workflow](https://github.com/XQP-Munich/LDPC4QKD/actions/workflows/ci-cmake_tests.yml/badge.svg)](https://github.com/XQP-Munich/LDPC4QKD/actions)
[![License](https://img.shields.io/github/license/XQP-Munich/LDPC4QKD)](./LICENSE)
# LDPC4QKD - LDPC codes for rate adaptive distributed source coding

Note: This repository is still missing some features which will be added very soon.

## Overview

This repository aims to solve the following problem. Suppose Alice and Bob each have a bit-string of length $N$. Bobs bit-string is the result of randomly flipping each bit of Alice with probability $p$, i.e., it is the output of a binary symmetric channel with **channel parameter** $p$. (Note: the same considerations apply for soft inputs instead of bit-strings and a Gaussian channel.) The goal is for Alice to send (via a noise-less channel) to Bob a message (called syndrome), such that Bob can recover Alice's bit-string given his bit-string and the syndrome received from Alice. There are two important metrics, which we want to minimize: 
- the probability that Bob will fail to recover Alice's bit-string, which is called the frame error rate (FER)
- the length of the syndrome

We define the syndrome to be the matrix-vector product (modulo 2) of a sparse binary matrix (LDPC matrix) and Alice's bit-string. In this case, Bob can use an inference algorithm (loopy belief propagation) to obtain a guess of Alice's bit-string. If the LDPC matrix has size $M \times N$, we call $M/N$ the **rate** of the LDPC matrix and $N$ the **block size**. (Note: the term rate is used differently in forward error correction, where it means $1 - M/N$.) Minimizing the length of the syndrome actually means minimizing the rate for a given channel parameter. Note that, for any given channel parameter, there is a trade-off between small rate, small FER and small block size. Furthermore, there are theoretical limits as to how small the rate can be (the Slepian Wolf limit, sometimes called Shannon limit, applies to assymptotically large bloc size. Limits for finite block sizes also exists but are more complicated).

### Contrast with forward error correction

Contrary to how forward error correction works, generator matrices for the LDPC code are not used at all for distributed source coding. This repository does not provide generator matrices (calculating them from the parity check matrices is straightforward).
Our decoder implementation operates on a bit-string (noisy version of true message) and its correct syndrome. This is slightly different from what is used for forward error correction (as in e.g. AFF3CT), where the decoder operates on only the noisy codeword. The noisy codeword is the result transmitting the codeword (true message encoded using a generator matrix) via a noisy channel.


## How to use

Start from the `examples` directory, which shows a basic example of how to use the C++ header containing the decoder. 
TODO

## How to contribute
This repository is actively maintained. Issues and pull requests will be responded to and processed.

Let us know if you're having problems with the provided materials, wish to contribute new ideas or need modifications or our work for your application.

## List of contents

### Data files
- A number of LDPC codes. Their parity check matrices are stored in a custom file format (called `CSCMAT`).
- Simulations results (done using AFF3CT) showing FER of the LDPC matrices at various channel parameters.
- Simulation results done using this repository only, showing FER of LDPC matrices, their rate adapted versions, and average rate under rate adaption (Work in progress!).
- For each LDPC matrix, a specification of rate adaption. This is a list of pairs of row indices of the matrix that are combined (added mod 2) in each rate adaption step.

### Julia code to
- load the LDPC matrices from `CSCMAT` files
- save the compressed sparse column (CSC) representation and also to export to other standard formats, such as `alist`.

### C++ code
- Basic LDPC decoder using belief propagation (BP). Contained in a single header file (`PATH`) and easy to use. Can perform syndrome computation as well. The LDPC code can be embedded into the executable or loaded from a file at program runtime.
- For applications that only require syndrome computation but no decoding, we provide a separate implementation for multiplication of a sparse binary matrix and a dense binary vector (LDPC syndrome computation). This is also contained in a single header file (`PATH`).
- For this purpose, the LDPC matrix can be stored within the executable. Julia code is provided to automatically generate header files storing an LDPC code in constant arrays.

TODO add paths

## Planned features and improvements
- Code to automatically generate reports on quality of all LDPC codes and rate adapted performance (work in progress)
- More and better LDPC codes with more sizes and rates
- The BP implementation is not state-of-the-art or heavily optimized. We may improve it in the future.
- Using compressed sparse row format may have advantages over the current compressed sparse column format. This needs to be tested.
- the quasi-cyclic structure of LDPC matrices can be exploited to store the matrices more efficiently. For the encoder, this may be combinable with storing individual bits as integers rather than booleans, again saving some memory. The encoding in this case would consist of adding bit-shifted integers.


## Appendix

### Table of LDPC matrices

See directory `codes`.
