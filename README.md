# LDPC4QKD

Rate adaptive LDPC-based constructions designed for distributed source coding.

Includes LDPC codes with specified rate adaption and example C++ implementation of encoder (just sparse matrix multiplication) and decoder (
Belief Propagation). Both decoder can be used as header-only C++ libraries.
The decoder can also encode; the dedicated encoder should only be used for applications where code size is critical.

Note: This repository is still missing several features, as well as the actual LDPC codes, which will be added soon.
