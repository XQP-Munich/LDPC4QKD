# This script generates a `.hpp` file (C++ header) containing
# an LDPC code stored in compressed sparse column (CSC) format.
# See command line help for how to use it by (in this directory) running
# `julia --project ldpc_codegen.jl --help`

using SparseArrays
using LinearAlgebra

# Must be installed (use the Pkg package manager to instantiate the proviced manifest)
using ArgParse
using LDPCStorage  # custom project. Helper methods for reading and writing LDPC matrices to files


struct ValidPath
    val::String
    ValidPath(s::AbstractString) = isfile(s) ? new(s) : error("No such file: '$s'")
end

"""
Parses command line arguments
"""
function parse_my_args()
    s = ArgParseSettings(
        "Create a C++ code-file containing an LDPC matrix in compressed sparse column (CSC) format.
        The matrix is given as an `.alist` file.
        Example usage:
        julia --project ldpc_codegen.jl --input_code_path ./filename.qccsc.json --output_path autogen_ldpc.hpp
        "
        ; allow_ambiguous_opts=false, preformatted_description=true)

    @add_arg_table! s begin
        "--input_code_path"
            required = true
            arg_type = ValidPath
            help = "Path to .alist or .json file containing the LDPC code.
                For .bincsc.json files, raw LDPC storage is supported.
                For .qccsc.json files, quasi-cyclic-exponents are supported"
        "--output_path"
            default = "autogen_ldpc_matrix_csc.hpp"
            arg_type = String
            help = "Path to put the automatically generated C++ file containing the LDPC code.
            If the file extension is .bincsc.json, a compressed sparse column storage of the full matrix
            (actual binary matrix, not QC-exponents) is written instead."
        "--output_quasi_cyclic_exponents"
            action = :store_true
            help = "ONLY SUPPORTED if `output_path` is a C++ header!!!
            If this flag is NOT set (default), the header contains two arrays (`colptr`, `row_idx`).
            If this flag is set, the resulting C++ header contains three arrays (`colptr`, `row_idx` and `values`).
            "
    end

    parsed_args = parse_args(s)
    println("Parsed args:")
    for (key, val) in parsed_args
        println("  $key  =>  $(repr(val))")
    end

    if parsed_args["output_quasi_cyclic_exponents"] && !endswith(abspath(parsed_args["output_path"]), ".hpp")
    @show parsed_input_code_abspath
        throw(ArgumentError(
                    "Flag `--output_quasi_cyclic_exponents` is only supported if output is a C++ header."))
    end
    parsed_input_code_abspath = abspath(parsed_args["input_code_path"].val)
    return parsed_input_code_abspath,
            abspath(parsed_args["output_path"]),
            parsed_args["output_quasi_cyclic_exponents"]
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


"""
Converts all `.qccsc.json` files in the folder given by `input_folder_path` to same name with extension `.bincsc.json`
files in the `output_folder_path`.
Each input file is assumed to contain QC exponents for an LDPC matrix (as defined by the `.qccsc.json` format).
The output files store the full binary LDPC matrix, which can be used as input to the C++ programs.

Warning: this may overwrite pre-existing files in the output folder.
"""
function convert_all_qc_to_binary(
    input_folder_path::AbstractString, output_folder_path::AbstractString;
    verbose=false
    )
    isdir(input_folder_path) || error("`$input_folder_path` is not a valid folder path.")
    isdir(output_folder_path) || mkdir(output_folder_path)

    file_paths = readdir(input_folder_path; join=true)
    filter!(p -> endswith(p, ".qccsc.json"), file_paths)
    verbose && @info "Converting $(length(file_paths)) files."

    for p in file_paths
        H = LDPCStorage.load_ldpc_from_json(p; expand_qc_exponents_to_binary=true)

        # store file with same name but replaced file extension reflecting contents
        out_file_path = joinpath(output_folder_path, basename(p)[1:end-11]*".bincsc.json")
        LDPCStorage.save_to_bincscjson(out_file_path, H)
        verbose && @info "Converted `$p` -> `$out_file_path`."
    end
    return nothing
end


function main()
    code_path, output_path, output_quasi_cyclic_exponents = parse_my_args()

    if endswith(code_path, ".alist")
        H = load_alist(code_path)
    elseif endswith(code_path, ".qccsc.json")
        # When outputting QC-exponents, do not expand the matrix to be binary!
        H = LDPCStorage.load_ldpc_from_json(code_path; expand_qc_exponents_to_binary=!output_quasi_cyclic_exponents)
    elseif endswith(code_path, ".bincsc.json")
        H = LDPCStorage.load_ldpc_from_json(code_path; expand_qc_exponents_to_binary=false)
    elseif endswith(code_path, ".cscmat")
        @warn "The `.cscmat` file format is deprecated. Consider converting to the json based formats instead!"
        H = load_cscmat_standard_or_qc_exponents(code_path)
    else
        throw(ArgumentError(
            "Unsupported input file at path $code_path.
            Expected file extension `qccsc.json`,`bincsc.json`, `.cscmat` (deprecated) or `.alist`."))
    end

    if endswith(output_path, ".alist")
        LDPCStorage.save_to_alist(output_path, H)
    elseif endswith(output_path, ".hpp")
        if output_quasi_cyclic_exponents
            expansion_factor = LDPCStorage.get_qc_expansion_factor(code_path)
            open(output_path, "w+") do io
                LDPCStorage.print_cpp_header_QC(io, H; expansion_factor)
            end
        else
            open(output_path, "w+") do io
                LDPCStorage.print_cpp_header(io, H)
            end
        end
    elseif endswith(output_path, ".cscmat")
        LDPCStorage.save_to_cscmat(H, output_path; allow_omit_entries_if_only_stored_ones=true,)
        # elseif code_path_file_extension == ".qccsc.json"  # need expansion factor. Could only do (useless) `qccsc.json` -> `qccsc.json`
    elseif endswith(output_path, ".bincsc.json")
        LDPCStorage.save_to_bincscjson(output_path, H)
    else
        throw(ArgumentError(
            "Unsupported output file format at path $output_path. Expected file extension `.alist`, `.hpp` `.bincsc.json` or `.cscmat`."))
    end

    @info "Saved output file ($(Base.stat(output_path).size / 1000.) kilobytes) at '$(output_path)'."
end


main()
