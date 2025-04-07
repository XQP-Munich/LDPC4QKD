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

The following is a list of LDPC matrices of size `rows x columns` that are provided in this repository.
The printed row and column counts are those of the actual binary parity check matrix.
The corresponding matrix of quasi-cyclic exponents, which is smaller by a factor of the quasi-cyclic expansion factor `QCE`.
We define the **rate** as the number of rows divided by the number of columns.

The tables list various file hashes (sha256) for the matrices.
These allow you to hash 
(using e.g. `sha256sum` [Linux](https://www.linux.org/docs/man1/sha256sum.html), or `Get-FileHash` [Windows PowerShell](https://learn.microsoft.com/en-us/powershell/module/microsoft.powershell.utility/get-filehash?view=powershell-7.4))
all files in the `/codes` directory and find the ones you need using the tables below.

- The sha256 hashes of **JSON files** connect the table items to JSON files in this repository.

- For frame error rate simulations, we use [AFF3CT](https://github.com/aff3ct/aff3ct), which accepts the `alist` format but not our custom format.
  To make sure which code is simulated, we include the sha256 of the **alist file** in the raw results (at `/codes/aff3ct_fer_simulations`)
  Alist files can be quite large (~50 Megabytes per file) and are therefore not part of this repository.
  To reproduce the FER results, use the command line interface or [LDPCStorage.jl](https://github.com/XQP-Munich/LDPCStorage.jl) to recreate alist files and check the hash.
  For plots and evaluation of the FER simulations, see the paper and Julia code at `/codes/aff3ct_fer_simulations/evaluation`.
  The Julia code generates all plots in the paper (and more) by reading all raw AFF3CT outputs and matching them by the sha256 of the respective alist file.

- The sha256 of corresponding **rate adaption files** (where applicable) connect each code to its rate adaption.

## Protograph-based quasi-cyclic matrices

| Rate |    rows / columns    |   Code Structure   |          Protograph          |          sha256 of alist file (NOT in this repository)           |             sha256 of json file (in this repository)             | sha256 of CORRESPONDING rate adaption file (in this repository)  |
|:----:|:--------------------:|:------------------:|:----------------------------:|:----------------------------------------------------------------:|:----------------------------------------------------------------:|:----------------------------------------------------------------:|
| 1/2  |    2,048 / 4,096     | Protograph QCE=32  |     `[1 2 1 3; 1 0 2 5]`     | 1fbeda66bd135033250aa88ef526f0bb5bb0a5dc9b61e7a960db1f03cb1dd935 | 098de6e117e43a408603758a3cb1985d9c18c188d08598485b22ab3b2235e8a5 | 21cb9c807dfe128035ca4784b6e737233e4b6ca2cc1f26ba5d815a2bb1f16e8e |
| 1/2  |    8,192 / 16,384    | Protograph QCE=64  |     `[1 2 1 3; 1 0 2 5]`     | 5bfa71c25ddb19f88a791fc15da9ecbe09dbe3bd49ebba87ecb596f5e1a6ea4f | 44dead953402ebe461f6c3895cc66b7f24366c6bd27ec84bb11de206778117a6 | e8b78b264d627da9ca4a53b0ffea4664a8ee2538a7b6a850c37799445a02b084 |
| 1/2  | 524,288 / 1,048,576  | Protograph QCE=512 |     `[1 2 1 3; 1 0 2 5]`     | 6f1747ed60f2956a03250282395baba2437d1684588cec7b58e63b395fe133ca | 9f8c301f67b663c673a6feec52c8cc8b122bef97fe9aa06208f634be2f652c6f | ab304590083fc586e0623a708ca612ccb80119803772beba26cd210ef5e3a756 |
| 1/3  |    2,048 / 6,144     | Protograph QCE=32  | `[3 1 3 4 2 2; 4 1 0 4 0 1]` | 54adc87fd548a4aa8c61efaf54194beca750afd72124ff52846bee4ee2cf482a | 608f5ab52838bf1c1660412824de51237879e2bc5a2b073369852a2d7c8a0c24 | 45e7db986efe7964d0e6f0c4fc55080eb663015a45dc63cb33780265827d9f2c |
| 1/3  |    8,192 / 24,576    | Protograph QCE=64  | `[3 1 3 4 2 2; 4 1 0 4 0 1]` | dd32e139f2ab999ec18d8c4933dcb112fbfa4a26b511f29f57cd71590c8440dc | 9b053763c2092802c02791da2d145b13bec8af11646d44bcb9f2db3284961606 | 3c84cc12a1be72f60335e5a131157814a132e91c84cb4a5fe19c0368d346e67a |
| 1/3  | 524,288  / 1,572,864 | Protograph QCE=512 | `[3 1 3 4 2 2; 4 1 0 4 0 1]` | d839b0af96478e8d1e6c80ce52236aa284fcffcdc6ef7ed1603598a5eb22f184 | 7f8df4cb9e4ef53d12f99813634e93e9d447ed7393d69e25dbd6e290ee601e43 | 2fbea43dad8dafac0363ab2c956c029b9b39626857334a6c31612bb3b4ec6d42 |

## Degree-distribution-based quasi-cyclic matrices

The degree distributions used to make an initial "proto-matrix" are from [arXiv:0901.2140v1](https://arxiv.org/abs/0901.2140).
The proto-matrix is then QC-lifted using the specified factor.
For the codes below, rate adaption is not stored in files but instead specified using an algorithm.

### TODO point towards explanation of rate adaption for these codes

| Rate | columns | Code Structure |             sha256 of json file (in this repository)             |          sha256 of alist file (NOT in this repository)           |
|------|:-------:|:--------------:|:----------------------------------------------------------------:|:----------------------------------------------------------------:|
| 0.50 | 819200  |    QCE=1024    | 9bc0e3aa3296969770eae484247c34154e731f9c722f3fb56c9a8a8fb7d65f19 | f4963997d779ad20a2aa15f924609d3b9bdca79ffc5cb8572f4594459ab7a4b4 |
| 0.45 | 819200  |    QCE=1024    | 0b4c996d4f75232494086756672f4f9cd4cd6529db04ce043f457be58b5b91b6 | f02a453b2152f445492b3130b5bc6c1ae3675710ff876eb6f67f08cfbe223af4 |
| 0.40 | 819200  |    QCE=1024    | 11841a98f507759ec918106931d8690bb4be16b7aaf9dfcec87a1a1eb567b896 | 34d30850b7a07170e9b1c2102445930c004f9dd1bcfaa383eb75873e46d355c5 |
| 0.35 | 819200  |    QCE=1024    | 1f9b7a72eba6225ba9991797cdfa015a0dbc2279f3299032e5d8e7221e1c199d | 18f69a856a9ab5ec94085039546624b161061f686b5aba01a51cbaf780434d8f |
| 0.30 | 819200  |    QCE=1024    | 5fa69f06625ab728b8902bb262d2f2a1f99c8a3104e8856af6bc48c419051f00 | 81e397215fefa7dbafd009ba419928730654562fc586c914b75c5fd4446b5a75 |
| 0.25 | 819200  |    QCE=1024    | d22d1cf8c91ae0eb454b4fe1e74dd156d58183bf0feb9800a50d0e53972aca64 | 03ea702290fd2e1e22d2d45cdf0400ef0ab90380543c5e5d63f9b7a0de1c0030 |
| 0.20 | 819200  |    QCE=1024    | f77ff0585710e66462a70b51748a6f0a19d4f3c2e90c181da88d91007ae6eb06 | 66100075e0ccfe345a235d1cc03e32bab7835e567a20ed212f5562088177ba44 |
| 0.15 | 819200  |    QCE=1024    | 8e55c523996bf20a7a6380d6f3730cc05886a4cd118c6613ae7bb160db8902fb | 82fca2da0043e0cf59b6cdf0bd73934cdb3ebcde64f5d11365a179056c1af245 |
| 0.10 | 819200  |    QCE=1024    | 5a4734315f1139259b3ef38f6f8d72d8533b848e97bc8c012143016e70e483d6 | d1f3a91a24bf72a2946de44b582bb846a74e180a22a6f38fe1ce99462a4a268f |

### Notes:

For protograph-based LDPC codes, the table gives the protograph as `[first row; second row; third row etc.]`. 

For quasi-cyclic (QC) codes, the quasi-cyclic expansion factor (or lifting degree) is given as QCE.
The QCE is the size of the block-sub-matrices constituting the full LDPC matrix defined by the quasi-cyclic exponents.
