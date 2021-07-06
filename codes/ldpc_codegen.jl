# This file generates a `.cpp` code file containing
# an LDPC code stored in compressed sparse column (CSC) format.
# See command line help for how to use it.

using SparseArrays
using LinearAlgebra

# Must be installed:
using ArgParse


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
        "alist_path"
            required = true
            arg_type = ValidPath
            help = "Path to .alist file containing LDPC code."
        "--output_path"
            default = "autogen_ldpc_matrix_csc.hpp"
            arg_type = String
            help = "Path to put the automatically generated C++ file containing the LDPC code."
        "--debug_use_only_part"
            help = "If true, only a very small part of the provided matrix is processed."
            action = :store_true
    end

    parsed_args = parse_args(args, s)
    println("Parsed args:")
    for (key, val) in parsed_args
        println("  $key  =>  $(repr(val))")
    end

    return abspath(parsed_args["alist_path"].val),
            abspath(parsed_args["output_path"]),
            parsed_args["debug_use_only_part"]
end


"""parse a single line of space separated integers"""
space_sep_ints(s::AbstractString) = parse.(Int, split(s))

 
function file_extension(path::String)
    if contains(path, ".")
        return path[findlast(isequal('.'),path):end]
    else
        return ""
    end
end 

 
 """Load an LDPC matrix from a text file in alist format."""
function load_alist(
    file_path::AbstractString; check_redundant=false
    )
    if file_extension(file_path) != ".alist"
        @warn "load_alist called on file with extension '$(file_extension(file_path))', expected '.alist'"
    end

    file = open(file_path, "r")
    nVN, nCN = space_sep_ints(readline(file))
    dmax_VN, dmax_CN = space_sep_ints(readline(file))
    var_node_degs = space_sep_ints(readline(file))
    check_node_degs = space_sep_ints(readline(file))
    remaining_lines = readlines(file)
    close(file)

    if length(remaining_lines) != nVN + nCN
        error("Number of lines in $file_path is inconcistent with stated matrix size.")
    end

    if dmax_CN != maximum(check_node_degs)
        error("Alist file $file_path claims: max. CN degree=$dmax_CN but contents give $(maximum(check_node_degs)).")
    end

    if dmax_VN != maximum(var_node_degs)
        error("Alist file $file_path claims: max. VN degree=$dmax_CN but contents give $(maximum(var_node_degs)).")
    end

    # parity check matrix
    H = spzeros(Int8, nCN, nVN)

    # fill the matrix
    for col_ind in 1:nVN
        rows = space_sep_ints(remaining_lines[col_ind])

        if check_redundant && length(rows) != var_node_degs[col_ind]
            error("Variable node degree in $file_path inconcistent with below data for VN $col_ind.")
        end

        for row_ind in rows
            H[row_ind, col_ind] = 1
        end
    end


    # the second half of the alist file is redundant. Check that it is consistent.
    if check_redundant
        entry_counter = 0
        for row_ind in 1:nCN
            cols = space_sep_ints(remaining_lines[nVN + row_ind])

            check_node_degree = length(cols)
            if check_node_degree != check_node_degs[row_ind]
                error("Check node degree in $file_path inconcistent with below data for CN $row_ind.")
            end

            entry_counter += check_node_degree
            for col_ind in cols
                if H[row_ind, col_ind] != 1
                    error("VN and CN specifications in $file_path disagree on matrix entry ($row_ind, $col_ind).")
                end
            end
        end

        if entry_counter != sum(H)
            error("VN and CN specification in $file_path are inconsistent.")
        end
    end

    return H
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

        println(f, "} // namespace RALDPC")
    end

end


function main(args)
    alist_path, output_path, only_debug_mode = parse_my_args(args)
    H = load_alist(alist_path)
    write_cpp_constexpr_CSC(H, output_path; only_debug_mode)

    @info "Saved C++ file ($(Base.stat(output_path).size) bytes) at '$(output_path)'."
end


main(ARGS)
