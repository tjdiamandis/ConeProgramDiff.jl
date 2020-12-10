
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
    file_vec = readlines("test.txt")
    nvars = length(file_vec)

    if dense
        ~(nvars in [3, 6]) && throw(ArgumentError("Invalid file"))
        offset = 0
        A_flat = parse.(Float64, split(file_vec[1], '\t'))
        A = reshape(A_flat, (m,n))
    else
        ~(nvars in [5, 8]) && throw(ArgumentError("Invalid file"))
        offset = 2
        A_row_inds = parse.(Int, split(file_vec[1], '\t'))
        A_col_inds = parse.(Int, split(file_vec[2], '\t'))
        A_vals = parse.(Float64, split(file_vec[3], '\t'))
        A = sparse(A_row_inds, A_col_inds, A_vals)
    end

    b = parse.(Float64, split(file_vec[2+offset], '\t'))
    c = parse.(Float64, split(file_vec[3+offset], '\t'))

    if nvars == 6
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
cp_from_file("test.txt", dense=false)
