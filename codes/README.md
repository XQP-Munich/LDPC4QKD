# Custom file formats
We use two custom formats to store LDPC codes as JSON files.

1. `.bincsc.json` files store a sparse binary LDPC matrix (matrix containing only ones and zeros) directly. 
The file contains two arrays (`rowidx` and `colptr`) of integers, which store the zero-based indices as used for compressed sparse column (CSC) storage of sparse matrices. 
The values array (used in CSC) is omitted due to all non-zero entries being ones.
For large matrices, this format leads to big files (about half the size of `.alist` files).
2. `.qccsc.json` files store a sparse matrix containing the exponents used to create a quasi-cyclic LDPC matrix.
This way of storing the matrix leads to much smaller files.
In this case, all three arrays defining the CSC storage are used. 
The values array (called `nzval`) stores the quasi-cyclic exponents.
The quasi-cyclic expansion factor is stored as `qc_expansion_factor`.
The numbers of rows and columns refer to the shape of the matrix of exponents.
Note that the C++ code does not yet support loading `.qccsc.json` directly.
Convert them to a c++ header file (for normal use) or to the `.bincsc.json` format (for simulations).

Our custom Julia package [LDPCStorage.jl](https://github.com/XQP-Munich/LDPCStorage.jl) contains functions to read and write such files.
They also support the `.alist` format and our deprecated `CSCMAT` format which was used by older versions of the project. 

# Storing LDPC matrix in C++ header
For use in the LDPC encoder/decoder, we store the LDPC matrices as static data in C++ header files and embed them into the executable.
We provide a command line interface to convert the JSON-based format to a C++ header file.
To do so, install Julia (at least v1.6) and perform these steps to download and install dependencies (as specified in `Project.toml` and `Manifest.toml`):

    cd LDPC4QKD/codes
    julia --project
    ] instantiate

Next, use these commands to convert LDPC codes to a C++ header file:

    cd LDPC4QKD/codes
    julia --project ldpc_codegen.jl --input_code_path <path>.json --output_path <path>.hpp


# List of LDPC matrices

The following is a list of LDPC matrices of size `rows x columns` that are provided in this repository. We define the **rate** as the numbers of rows divided by the number of columns.

The sha256 hashes of json files connect the table items to files in this reporitory. The sha256 hashes of alist files are included in raw [AFF3CT](https://github.com/aff3ct/aff3ct) frame error rate (FER) simulation outputs (folder `codes/aff3ct_fer_simulations`, does not use rate adaption) to make sure they refer to the same code. 
([AFF3CT](https://github.com/aff3ct/aff3ct) accepts `alist` format but not our custom format. 
Alist files can be quite large and are therefore not part of this repository. 
To reproduce the [AFF3CT](https://github.com/aff3ct/aff3ct) FER results, use the command line interface or [LDPCStorage.jl](https://github.com/XQP-Munich/LDPCStorage.jl) to create alist files.)

| Rate |    rows / columns    |   Code Structure   | Protograph or Degree Distribution |                       sha256 of alist file                       |                       sha256 of json file                        |
|:----:|:--------------------:|:------------------:|:---------------------------------:|:----------------------------------------------------------------:|:----------------------------------------------------------------:|
| 1/2  |    2,048 / 4,096     | Protograph QCE=32  |       `[1 2 1 3; 1 0 2 5]`        | 1fbeda66bd135033250aa88ef526f0bb5bb0a5dc9b61e7a960db1f03cb1dd935 | 098de6e117e43a408603758a3cb1985d9c18c188d08598485b22ab3b2235e8a5 |
| 1/2  |    8,192 / 16,384    | Protograph QCE=64  |       `[1 2 1 3; 1 0 2 5]`        | 5bfa71c25ddb19f88a791fc15da9ecbe09dbe3bd49ebba87ecb596f5e1a6ea4f | 44dead953402ebe461f6c3895cc66b7f24366c6bd27ec84bb11de206778117a6 |
| 1/2  | 524,288 / 1,048,576  | Protograph QCE=512 |       `[1 2 1 3; 1 0 2 5]`        | 6f1747ed60f2956a03250282395baba2437d1684588cec7b58e63b395fe133ca | 9f8c301f67b663c673a6feec52c8cc8b122bef97fe9aa06208f634be2f652c6f |
| 1/3  |    2,048 / 6,144     | Protograph QCE=32  |   `[3 1 3 4 2 2; 4 1 0 4 0 1]`    | 54adc87fd548a4aa8c61efaf54194beca750afd72124ff52846bee4ee2cf482a | 608f5ab52838bf1c1660412824de51237879e2bc5a2b073369852a2d7c8a0c24 |
| 1/3  |    8,192 / 24,576    | Protograph QCE=64  |   `[3 1 3 4 2 2; 4 1 0 4 0 1]`    | dd32e139f2ab999ec18d8c4933dcb112fbfa4a26b511f29f57cd71590c8440dc | 9b053763c2092802c02791da2d145b13bec8af11646d44bcb9f2db3284961606 |
| 1/3  | 524,288  / 1,572,864 | Protograph QCE=512 |   `[3 1 3 4 2 2; 4 1 0 4 0 1]`    | d839b0af96478e8d1e6c80ce52236aa284fcffcdc6ef7ed1603598a5eb22f184 | 7f8df4cb9e4ef53d12f99813634e93e9d447ed7393d69e25dbd6e290ee601e43 |

### Notes:

For Protograph-based LDPC codes, the table gives the protograph as `[first row; second row; third row etc.]`. 

For quasi-cyclic (QC) codes, the quasi-cyclic expansion factor (or lifting degree) is given as QCE.
The QCE is the size of the block-sub-matrices constituting the full LDPC matrix defined by the quasi-cyclic exponents.
