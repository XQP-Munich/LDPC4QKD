# LDPC matrices

The following is a list of LDPC matrices of size MxN that are provided in this repository. The sha256 hashes of data files are given to make clear which data files and simulations refer to which code. The alist file hash is given to connect it to AFF3CT simulation results (AFF3CT accepts alist format but not our custom format. Alist files can be quite large and are therefore not part of this repository. To reproduce the aff3ct simulation results, Julia code is provided to convert from cscmat to alist format.)

| Rate | M (= Rate * N) |    N    | Code Structure |    sha256 of alist file     |    sha256 of cscmat file (exponents)   |
|------|----------------|---------|----------------|----------------------------|-----------------------------|
| 1/2  | 2048           |   4096 | Protograph-QC  | 1fbeda66bd135033250aa88ef526f0bb5bb0a5dc9b61e7a960db1f03cb1dd935                        | 12cdb1acbe918b2db8efce2c897dcd0ccb3ae9a4af98220713f199eec0c874d3                             |
| 1/2  | 8192           | 16384 | Protograph-QC  | 5bfa71c25ddb19f88a791fc15da9ecbe09dbe3bd49ebba87ecb596f5e1a6ea4f                        | 2207bee57d8c8e05fabdeea6585e476f0dcbfa89f37fcfed1374b9ade13dbe12  |
| 1/2  | 524288         | 1048576 | Protograph-QC  | 6f1747ed60f2956a03250282395baba2437d1684588cec7b58e63b395fe133ca                        | 98e9fc7b26822043c894ad6c842e823278c317c958ffafe17179bc0124f85ee7                             |
| 1/3  | 2048           | 6144 | Protograph-QC  | 54adc87fd548a4aa8c61efaf54194beca750afd72124ff52846bee4ee2cf482a                        | f40c5d91891e54f5ad44d58fc8fb970bd379af829cb6c9eb58eb546f00c6c91b                             |
| 1/3  | 8192           | 24576 | Protograph-QC  | d839b0af96478e8d1e6c80ce52236aa284fcffcdc6ef7ed1603598a5eb22f184                        | 5502076bac2654824b58fe1744d106341b97c4f0c03c1be001d2f9bff07f273b                             |
| 1/3  | 524288         | 1572864 | Protograph-QC  | dd32e139f2ab999ec18d8c4933dcb112fbfa4a26b511f29f57cd71590c8440dc                        | ba59da531aa7683ee0a6ccc913d2dc58b449c6b0345acdb565ff2fc1bbfac962                             |

