module StorageUtilsLDPC

using SparseArrays
using LinearAlgebra

export load_alist, save_to_alist, save_to_cscmat, load_cscmat


"""Load an LDPC matrix from a text file in alist format."""
function load_alist(file_path::AbstractString; check_redundant=false,)
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


"""
    function save_to_alist(matrix::AbstractArray{Int8,2}, out_file_path::String)

Save LDPC matrix to file in alist format. For details about the format, see:
https://aff3ct.readthedocs.io/en/latest/user/simulation/parameters/codec/ldpc/decoder.html#dec-h-path-image-required-argument
http://www.inference.org.uk/mackay/codes/alist.html
todo test this carefully
"""
function save_to_alist(matrix::AbstractArray{Int8,2}, out_file_path::String)

    (the_M, the_N) = size(matrix)

    variable_node_degrees = get_variable_node_degrees(matrix)
    check_node_degrees = get_check_node_degrees(matrix)

    # write data as specified by the alist format
    lines = String[]
    # -- Part 1 --
    # 'the_N' is the total number of variable nodes and 'the_M' is the total number of check nodes
    push!(lines, "$the_N $the_M")

    # 'dmax_VN' is the highest variable node degree and 'dmax_CN' is the highest check node degree
    push!(lines, "$(maximum(variable_node_degrees)) $(maximum(check_node_degrees))")

    # list of the degrees for each variable node
    push!(lines, join(["$deg" for deg in variable_node_degrees], " "))

    # list of the degrees for each check node
    push!(lines, join(["$deg" for deg in check_node_degrees], " "))

    # -- Part 2 --
    # each following line describes the check nodes connected to a variable node, the first
    # check node index is '1' (and not '0')
    # variable node '1'
    """
        Get indices of elements equal to one in a matrix.
        :param matrix: 2d numpy array
        :return: Array of strings, one string with indices for each line in the matrix.
    """
    function get_node_indices(matrix::AbstractArray{Int8,2})
        # the first check node index is '1' (and not '0')
        degrees = [findall(row .== 1) for row in eachrow(matrix)]
        return [join(string.(index_list), " ") for index_list in degrees]
    end
    append!(lines, get_node_indices(transpose(matrix)))

    # -- Part 3 --
    # each following line describes the variables nodes connected to a check node, the first
    # variable node index is '1' (and not '0')
    # check node '1'
    append!(lines, get_node_indices(matrix))

    open(out_file_path, "w+") do file
        for line in lines
            println(file, line)
        end
    end

    return nothing
end



 # helper methods for testing the parity check matrix
 function get_variable_node_degrees(matrix::AbstractArray{Int8,2})
    @assert(length(size(matrix)) == 2, "Matrix required. Wrong number of dimensions")
    return [sum(row) for row in eachcol(matrix)]
 end


 function get_check_node_degrees(matrix::AbstractArray{Int8,2})
    @assert(length(size(matrix)) == 2, "Matrix required. Wrong number of dimensions")
    return [sum(row) for row in eachrow(matrix)]
 end


# Methods relevant to CSCMAT format

CSCMAT_FORMAT_VERSION = v"0.0.2"  # track version of our custom CSCMAT file format.


"""
write the three arrays defining compressed sparse column (CSC) storage of a matrix into a file.

If `try_hex`, integers in arrays are stored as hexadecimals (without 0x prefix!)
If `allow_omit_entries_if_only_stored_ones`, the `stored values` array is omitted if all stored values compare equal to 1.
"""
function save_to_cscmat(
    mat::SparseMatrixCSC, destination_file_path::String
    ;
    additional_header_lines="",
    try_hex::Bool=false,
    allow_omit_entries_if_only_stored_ones=false,
    )

    try
        additional_header_lines = "#"*join(split(additional_header_lines, "\n"), "\n# ")
    catch
        @warn "Failed to process additional header lines. Discarding them."
        additional_header_lines = ""
    end

    if file_extension(destination_file_path) != ".cscmat"
        @warn "Writing to sparse column storage file with with extension
            '$(file_extension(destination_file_path))', expected '.cscmat'. (path '$(destination_file_path)')"
    end

    number_map(x::Real) = x
    number_map(n::Integer) = string(n, base=try_hex ? 16 : 10)

    open(destination_file_path, "w+") do file
        println(file, "# $CSCMAT_FORMAT_VERSION")
        println(file, "# Compressed sparse column storage of matrix (arrays `colptr`, `rowval`, `stored_values`"
            *" as space separated $(try_hex ? "hexadecimal" : "decimal") integers. Stored entries may be zero.).")
        println(file, additional_header_lines)
        println(file, "# n_rows n_columns n_stored_entries")
        println(file, "$(mat.m) $(mat.n) $(nnz(mat))\n")

        println(file, join(number_map.(mat.colptr .- 1), " "))  # convert to zero-based indices
        println(file, "")
        println(file, join(number_map.(mat.rowval .- 1), " "))  # convert to zero-based indices
        println(file, "")

        if !allow_omit_entries_if_only_stored_ones || any(x->x!=1, mat.nzval)
            println(file, join(number_map.(mat.nzval), " "))
        end
    end

    return nothing
end


"""
read the three arrays defining compressed sparse column (CSC) storage of a matrix into a file.

If `try_hex`, integers are stored as hexadecimals (without 0x prefix!)
"""
function load_cscmat(file_path::String;
    print_file_header=false)
    expected_file_extension = ".cscmat"
    if file_extension(file_path) != expected_file_extension
        @warn "load_cscmat called on file '$(file_path)'
            with extension '$(file_extension(file_path))', expected $expected_file_extension"
    end

    file = open(file_path, "r")
    header = ""
    next_line = ""
    while true
        header *= (next_line*"\n")
        next_line = readline(file)

        (length(next_line) > 0 && next_line[1] == '#') || break
    end

    try
        header[1] == "# $CSCMAT_FORMAT_VERSION" && @warn "File written in format $(header[1][2:end]) is being read in format $CSCMAT_FORMAT_VERSION"
    catch e
        @warn "Failed to verify CSC file format version: $e"
    end

    if contains(header, "hexadecimal")
        base = 16
    else
        base = 10
    end

    print_file_header && print(header)

    n_rows, n_cols, n_nnz = space_sep_ints(next_line)
    _ = readline(file)  # empty line

    colptr = space_sep_ints(readline(file); base)
    _ = readline(file)  # empty line

    rowval = space_sep_ints(readline(file); base)
    _ = readline(file)  # empty line

    stored_entries = space_sep_ints(readline(file); base)
    _ = readline(file)  # empty line
    remaining_lines = readlines(file)
    close(file)

    if length(remaining_lines) > 0
        @warn "Ignoring additional lines:\n`$remaining_lines`"
    end

    if length(stored_entries) == 0
        stored_entries = ones(Int8, n_nnz)
    end

    # convert from zero-based to one-based indices.
    return SparseMatrixCSC(n_rows, n_cols, colptr .+1, rowval .+1, stored_entries)
end


"""parse a single line of space separated integers"""
space_sep_ints(s::AbstractString; base=10) = parse.(Int, split(s); base)


function file_extension(path::String)
    if contains(path, ".")
        return path[findlast(isequal('.'),path):end]
    else
        return ""
    end
end


end  # module StorageUtilsLDPC
