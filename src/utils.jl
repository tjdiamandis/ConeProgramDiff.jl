
# Vector <-> matrix
vec_symm = MOSD.vec_symm
unvec_symm = MOSD.unvec_symm


# Write and read CP from file
function cp_to_file(file, params; opt_vals=(), dense=true)
    if dense
        writedlm(file, (params..., opt_vals...))
    else
        writedlm(file, (findnz(params[1])..., params[2:3]..., opt_vals...))
    end
end

function cp_from_file(file; dense=true)
    file_vec = readlines(file)
    nvars = length(file_vec)

    if dense
        ~(nvars in [3, 6]) && throw(ArgumentError("Invalid file"))
        offset = 0
        A_flat = parse.(Float64, split(file_vec[1], '\t'))
        b = parse.(Float64, split(file_vec[2+offset], '\t'))
        c = parse.(Float64, split(file_vec[3+offset], '\t'))
        m,n = length(b), length(c)
        A = reshape(A_flat, (m,n))
    else
        ~(nvars in [6, 8]) && throw(ArgumentError("Invalid file"))
        offset = 2
        A_row_inds = parse.(Int, split(file_vec[1], '\t'))
        A_col_inds = parse.(Int, split(file_vec[2], '\t'))
        A_vals = parse.(Float64, split(file_vec[3], '\t'))
        A = sparse(A_row_inds, A_col_inds, A_vals)
        b = parse.(Float64, split(file_vec[2+offset], '\t'))
        c = parse.(Float64, split(file_vec[3+offset], '\t'))
    end

    if nvars in [6, 8]
        x_star = parse.(Float64, split(file_vec[4+offset], '\t'))
        y_star = parse.(Float64, split(file_vec[5+offset], '\t'))
        s_star = parse.(Float64, split(file_vec[6+offset], '\t'))
    else
        x_star, y_star, s_star = nothing, nothing, nothing
    end
    return Dict(:A => A, :b => b, :c => c, :x_star => x_star,
                :y_star => y_star, :s_star => s_star)
end


# A = sparse(randn(4,3))
# b = randn(4)
# c = randn(3)
# cp_to_file("test.txt", (A, b, c), dense=false)
folder = "/home/csquires/Desktop"
program_sparse = cp_from_file("$folder/test_sparse_py.txt", dense=false)
program_dense = cp_from_file("$folder/test_dense_py.txt", dense=true)
sparse_params = (program_sparse[:A], program_sparse[:b], program_sparse[:c])
sparse_optvals = (program_sparse[:x_star], program_sparse[:y_star], program_sparse[:s_star])
cp_to_file("$folder/test_sparse_jl.txt", sparse_params, opt_vals=sparse_optvals, dense=false)
dense_params = (program_dense[:A], program_dense[:b], program_dense[:c])
dense_optvals = (program_dense[:x_star], program_dense[:y_star], program_dense[:s_star])
cp_to_file("$folder/test_dense_jl.txt", dense_params, opt_vals=dense_optvals, dense=true)
