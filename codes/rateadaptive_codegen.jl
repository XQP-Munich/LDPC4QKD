# This script generates a `.hpp` file (C++ header) containing
# an LDPC code stored in compressed sparse column (CSC) format.
# See command line help for how to use it.

using SparseArrays
using LinearAlgebra
using DelimitedFiles

using ArgParse # Must be installed (use the Pkg package manager)


header = raw"""
// This is an automatically generated file.
// It contains row indices that should be combined to reduce the number of syndrome bits
// of an LDPC code in a rate-adaptive manner.
// These row indices are obtained via optimization of the particular LDPC matrix.

#include <cstdint>
#include <array>

namespace AutogenRateAdapt {

"""


struct ValidPath
    val::String
    ValidPath(s::AbstractString) = isfile(s) ? new(s) : error("No such file: '$s'")
end


function parse_my_args(args)
    s = ArgParseSettings(
        "Create a C++ code-file containing line indices for combination, 
        to be applied to a given LDPC code in order to achieve rate adaption."
        ; allow_ambiguous_opts=false)

    @add_arg_table! s begin
        "rate_adaption_file_path"
            required = true
            arg_type = ValidPath
            help = "Path to file containing the line indices for rate adaption."
        "--output_path"
            default = "autogen_rate_adaption.hpp"
            arg_type = String
            help = "Path to put the automatically generated C++ file containing the rate adaption."
        "--debug_use_only_part"
            help = "If true, only a small part of the provided data processed."
            action = :store_true
    end

    parsed_args = parse_args(args, s)
    println("Parsed args:")
    for (key,val) in parsed_args
        println("  $key  =>  $(repr(val))")
    end

    return abspath(parsed_args["rate_adaption_file_path"].val),
            abspath(parsed_args["output_path"]),
            parsed_args["debug_use_only_part"]
end


function write_cpp_constexpr_rate_adaption(
    rate_adaption_rows::AbstractArray{T, 2} where T <: Integer,
    output_cpp_path::AbstractString
    ; only_debug_mode = false
    )
    size(rate_adaption_rows, 2) == 2 || throw(
        ArgumentError("Expecting input of shape N x 2, where N is the number of time that two rows are combined.
                        Received shape $(size(rate_adaption_rows))"))

    if only_debug_mode  # make input data smaller...
        cut_at = minimum([10, size(rate_adaption_rows)...])
        rate_adaption_rows = rate_adaption_rows[1:cut_at, :]
    end

    max_row_idx = maximum(rate_adaption_rows)
    if log2(max_row_idx) < 16
        row_idx_type = "std::uint16_t"
    elseif log2(max_row_idx) < 32
        row_idx_type = "std::uint32_t"
    elseif log2(max_row_idx) < 64
        row_idx_type = "std::uint64_t"
    else
        throw(ArgumentError("Input insanely large? Has $max_row_idx entries..."))
    end

    open(output_cpp_path, "w+") do f
        print(f, header)

        println(f, """
        constexpr std::size_t num_combined_rows = $(length(rate_adaption_rows));\n
        constexpr std::array<$row_idx_type, num_combined_rows> rows = {""")

        for (i, idx) in enumerate(transpose(rate_adaption_rows))
            print(f, "0x$(string(idx - 1, base=16))")  # Convert index to base zero
            if i != length(rate_adaption_rows)
                print(f, ",")
            end
            if mod(i, 100) == 0
                println(f, "")  # for formatting.
            end
        end
        println(f, "\n};\n")

        println(f, "} // namespace RALDPC")
    end
end


function read_rate_adaption_file(file_path::AbstractString)
    return DelimitedFiles.readdlm(file_path, ',', Int)
end


function write_array_to_csv(arr::AbstractArray, out_path::AbstractString)
    open(out_path, "w") do io
       writedlm(io, arr, ',')
    end
end


function main(args)
    rate_adaption_file_path, output_path, only_debug_mode = parse_my_args(args)
    rate_adaption_rows = read_rate_adaption_file(rate_adaption_file_path)
    write_cpp_constexpr_rate_adaption(rate_adaption_rows, output_path; only_debug_mode)

    @info "Saved C++ file ($(Base.stat(output_path).size) bytes) at '$(output_path)'."
end


main(ARGS)
