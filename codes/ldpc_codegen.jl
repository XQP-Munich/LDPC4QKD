# This script generates a `.hpp` file (C++ header) containing
# an LDPC code stored in compressed sparse column (CSC) format.
# See command line help for how to use it.

using SparseArrays
using LinearAlgebra

# Must be installed (use the Pkg package manager to instantiate the proviced manifest)
using ArgParse
using LDPCStorage  # custom project. Helper methods for reading and writing LDPC matrices to files


header = raw"""
// This is an automatically generated file.
// A sparse LDPC matrix (containing only zeros and ones) is saved in compressed sparse column (CSC) format.
// Since the matrix (and LDPC code) is known at compile time, there is no need to save it separately in a file.
// This significantly blows up the executable size (the memory would still have to be used when saving the matrix).
// The method seems to be reasonably fast (on a standard laptop).

#include <cstdint>
#include <array>

namespace AutogenLDPC {

"""


struct ValidPath
    val::String
    ValidPath(s::AbstractString) = isfile(s) ? new(s) : error("No such file: '$s'")
end


function parse_my_args(args)
    s = ArgParseSettings(
        "Create a C++ code-file containing an LDPC matrix in compressed sparse column (CSC) format.
        The matrix is given as an `.alist` file."
        ; allow_ambiguous_opts=false)

    @add_arg_table! s begin
        "--input_code_path"
            required = true
            arg_type = ValidPath
            help = "Path to .alist or .cscmat file containing the LDPC code. 
                For .cscmat files, both qc-exponents as well as raw LDPC storage is supported."
        "--output_path"
            default = "autogen_ldpc_matrix_csc.hpp"
            arg_type = String
            help = "Path to put the automatically generated C++ file containing the LDPC code.
            If the file extension is .cscmat, a compressed sparse column storage of the full matrix
            (actual binary matrix, not QC-exponents) is written instead."
        "--debug_use_only_part"
            help = "If true, only a very small part of the provided matrix is processed."
            action = :store_true
    end

    parsed_args = parse_args(args, s)
    println("Parsed args:")
    for (key, val) in parsed_args
        println("  $key  =>  $(repr(val))")
    end

    return abspath(parsed_args["input_code_path"].val),
            abspath(parsed_args["output_path"]),
            parsed_args["debug_use_only_part"]
end


function write_cpp_constexpr_CSC(
    H::AbstractArray{Int8, 2}, 
    output_cpp_path::AbstractString
    ; only_debug_mode = false
    )
    if only_debug_mode  # make matrix smaller...
        cut_at = minimum([10, size(H)...])
        H = sparse(H[1:ceil(Int, cut_at/2), 1:cut_at])
    else
        H = sparse(H)
    end

    _, _, values = findnz(H)

    all(values .== 1) || throw(ArgumentError("Expected matrix containing only zeros and ones."))

    num_nonzero = length(values)
    if log2(num_nonzero) < 16
        colptr_cpp_type = "std::uint16_t"
    elseif log2(num_nonzero) < 32
        colptr_cpp_type = "std::uint32_t"
    elseif log2(num_nonzero) < 64
        colptr_cpp_type = "std::uint64_t"
    else
        throw(ArgumentError("Input matrix not sparse? Has $num_nonzero entries..."))
    end

    if log2(size(H, 1)) < 16
        row_idx_type = "std::uint16_t"
    else 
        row_idx_type = "std::uint32_t"
    end

    open(output_cpp_path, "w+") do f
        print(f, header)

        N = size(H, 2)

        println(f, """
        constexpr std::size_t M = $(size(H, 1));
        constexpr std::size_t N = $(size(H, 2));
        constexpr std::size_t num_nz = $num_nonzero;

        constexpr std::array<$colptr_cpp_type, N + 1> colptr = {""")

        for (i, idx) in enumerate(H.colptr)
            print(f, "0x$(string(idx - 1, base=16))")  # Convert index to base zero
            if i != length(H.colptr)
                print(f, ",")
            end
            if mod(i, 100) == 0
                println(f, "")  # for formatting.
            end
        end
        println(f, "\n};\n")

        println(f, "// ------------------------------------------------------- \n")
        println(f, "constexpr std::array<$row_idx_type, num_nz> row_idx = {")

        for (i, idx) in enumerate(H.rowval)
            print(f, "0x$(string(idx - 1, base=16))")  # Convert index to base zero
            if i != length(H.rowval)
                print(f, ",")
            end
            if mod(i, 100) == 0
                println(f, "")  # for formatting.
            end
        end
        println(f, "\n};\n\n")

        println(f, "} // namespace AutogenLDPC")
    end

end


function load_cscmat_standard_or_qc_exponents(code_path; mode=:AUTO)
    if mode == :AUTO
        header = LDPCStorage.read_file_header(code_path)
        if contains(header, "Quasi cyclic exponents")
            mode = :QC_EXPONENTS
        else
            mode = :UNSTRUCTURED
        end
    end

    if mode == :QC_EXPONENTS
        return load_matrix_from_qc_cscmat_file(code_path)
    elseif mode == :UNSTRUCTURED
        H = load_cscmat(code_path)
        if any(x-> 0 <= x <= 1, H)
            @warn "Failed to obtain binary LDPC matrix from CSCMAT file!"
        end
        return H
    else
        @error "unsupported mode $mode"
    end
end


function main(args)
    code_path, output_path, only_debug_mode = parse_my_args(args)
    if LDPCStorage.file_extension(code_path) == ".alist"
        H = load_alist(code_path)
    elseif LDPCStorage.file_extension(code_path) == ".cscmat"

        H = load_cscmat_standard_or_qc_exponents(code_path)
    else
        @error "Unsupported input file at path $code_path. Expected file extension `.cscmat` or `.alist`."
    end

    if LDPCStorage.file_extension(output_path) == ".hpp"
        write_cpp_constexpr_CSC(H, output_path; only_debug_mode)
    elseif LDPCStorage.file_extension(output_path) == ".cscmat"
        save_to_cscmat(H, output_path; allow_omit_entries_if_only_stored_ones=true,)
    else
        @error "Unsupported input file at path $output_path. Expected file extension `.hpp` or `.cscmat`."
    end


    @info "Saved output file ($(Base.stat(output_path).size) bytes) at '$(output_path)'."
end


main(ARGS)
