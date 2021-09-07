# LDPC matrices and Julia code

## Julia code to load and pre-process LDPC codes
Use the Julia code contained in this directory to load and manipulate the LDPC codes. Functions are provided to load and save LDPC matrices in `.alist` and `.cscmat` file formats.

## CSCMAT files
We use a custom file format (file extension `.cscmat`) to store LDPC codes. Generally, such a file can contain any sparse matrix. We use the format for two distinct scenarios:

1. `QC-Format`: stores a sparse binary LDPC matrix (matrix containing only ones and zeros) directly. In this case, the file contains two arrays of integers, which store the zero-based indices of `rowidx` and `colptr`, as commonly used for compressed sparse column (CSC) storage of sparse matrices. The values array is omitted due to non-zero entries all being ones.
2. `Binary-Format`: stores a sparse matrix containing the exponents used to create a quasi-cyclic LDPC matrix. In this case, all three arrays defining the CSC storage are used. The `values` array stores the quasi-cyclic exponents.
   
Which of the two formats a given `.cscmat` file follows is clear from both the number of arrays (lines of space separated integers) contained in the text-file, as well as its header. Format 2 is used by default for the LDPC matrices provided in this folder (`codes`). This is because it saves a lot of disk space and allows storing a lot of very large matrices.

Use the Julia code to convert the array of quasi-cyclic components into either 
1. a C++ header file (storing the `rowidx` and `colptr` of the full binary LDPC matrix as constexpr arrays)
2. A `.cscmat` file in the `Binary-Format`. Such files can be imported by the C++ code from the file at program runtime and be used to initialize the decoder.

## Converting CSCMAT files

The Julia code provided in this repository can be called with command line options. To do so, first install Julia v1.6 or later. Then, perform these steps to install dependencies:

    cd LDPC4QKD/codes
    julia --project
    ] instantiate

After the installation is complete, use these commands to convert LDPC codes:

1. Convert to a C++ header file (default name `autogen_ldpc_matrix_csc.hpp`)

        cd LDPC4QKD/codes
        julia --project ldpc_codegen.jl --input_code_path <path>.cscmat --output_path <path>.hpp

2. Convert to a CSCMAT file that can be read by the C++ decoder

        cd LDPC4QKD/codes
        julia --project ldpc_codegen.jl --input_code_path <path>.cscmat --output_path <path>.cscmat

# List of LDPC matrices

The following is a list of LDPC matrices of size MxN that are provided in this repository. The sha256 hashes of data files are given to make clear which data files and simulations refer to which code. The alist file hash is given to connect it to [AFF3CT](https://github.com/aff3ct/aff3ct) simulation results ([AFF3CT](https://github.com/aff3ct/aff3ct) accepts alist format but not our custom format. Alist files can be quite large and are therefore not part of this repository. To reproduce the [AFF3CT](https://github.com/aff3ct/aff3ct) simulation results, Julia code is provided to convert from cscmat to alist format.)

| Rate | M (= Rate * N) |    N    | Code Structure |    sha256 of alist file     |    sha256 of cscmat file (exponents)   |
|------|----------------|---------|----------------|----------------------------|-----------------------------|
| 1/2  | 2048           |   4096 | Protograph-QC  | 1fbeda66bd135033250aa88ef526f0bb5bb0a5dc9b61e7a960db1f03cb1dd935                        | 12cdb1acbe918b2db8efce2c897dcd0ccb3ae9a4af98220713f199eec0c874d3                             |
| 1/2  | 8192           | 16384 | Protograph-QC  | 5bfa71c25ddb19f88a791fc15da9ecbe09dbe3bd49ebba87ecb596f5e1a6ea4f                        | 2207bee57d8c8e05fabdeea6585e476f0dcbfa89f37fcfed1374b9ade13dbe12  |
| 1/2  | 524288         | 1048576 | Protograph-QC  | 6f1747ed60f2956a03250282395baba2437d1684588cec7b58e63b395fe133ca                        | 98e9fc7b26822043c894ad6c842e823278c317c958ffafe17179bc0124f85ee7                             |
| 1/3  | 2048           | 6144 | Protograph-QC  | 54adc87fd548a4aa8c61efaf54194beca750afd72124ff52846bee4ee2cf482a                        | f40c5d91891e54f5ad44d58fc8fb970bd379af829cb6c9eb58eb546f00c6c91b                             |
| 1/3  | 8192           | 24576 | Protograph-QC  | d839b0af96478e8d1e6c80ce52236aa284fcffcdc6ef7ed1603598a5eb22f184                        | 5502076bac2654824b58fe1744d106341b97c4f0c03c1be001d2f9bff07f273b                             |
| 1/3  | 524288         | 1572864 | Protograph-QC  | dd32e139f2ab999ec18d8c4933dcb112fbfa4a26b511f29f57cd71590c8440dc                        | ba59da531aa7683ee0a6ccc913d2dc58b449c6b0345acdb565ff2fc1bbfac962                             |

[AFF3CT](https://github.com/aff3ct/aff3ct) simulations for all LDPC codes (without considering rate adaption) can be found in the folder `codes/aff3ct_fer_simulations` and are marked with the sha256 hash of the corresponding alist file.
