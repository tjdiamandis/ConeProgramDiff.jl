# https://github.com/cvxgrp/scs
"""
   unvec_symm(x, dim)
Returns a dim-by-dim symmetric matrix corresponding to `x`.
`x` is a vector of length dim*(dim + 1)/2, corresponding to a symmetric matrix
```
X = [ X11     X12/√2 ... X1k/√2
      X21/√2  X22    ... X2k/√2
      ...
      Xk1/√2  Xk2/√2 ... Xkk ],
```
where
`vec(X) = (X11, X12, X22, X13, X23, ..., Xkk)`

Note that the factor √2 preserves inner products:
`x'*c = Tr(unvec_symm(c, dim) * unvec_symm(x, dim))`
"""
function unvec_symm_scs(x, dim)
    X = zeros(eltype(x), dim, dim)
    idx = 1
    for i in 1:dim
        for j in 1:i
            if i == j
                X[i,j] = x[idx]
            else
                X[j,i] = X[i,j] = x[idx] / sqrt(2)
            end
            idx += 1
        end
    end
    return X
end


"""
   vec_symm(X)
Returns a vectorized representation of a symmetric matrix `X`.
`vec(X) = (X11, √2*X12, X22, √2*X13, X23, ..., Xkk)`

Note that the factor √2 preserves inner products:
`x'*c = Tr(unvec_symm(c, dim) * unvec_symm(x, dim))`
"""
function vec_symm_scs(X)
    x_vec = sqrt(2).*X[LinearAlgebra.tril(trues(size(X)))']
    idx = 1
    for i in 1:size(X)[1]
        x_vec[idx] =  x_vec[idx]/sqrt(2)
        idx += i + 1
    end
    return x_vec
end


# Write and read CP from file
function cp_to_file(file, params_dict; dense=true)
    A, b, c = params_dict[:A], params_dict[:b], params_dict[:c]
    if isnothing(params_dict[:x_star])
         _cp_to_file(file, (A, b, c), dense=dense)
    else
        x_star = params_dict[:x_star]
        y_star = params_dict[:y_star]
        s_star = params_dict[:s_star]
        _cp_to_file(file, (A, b, c), opt_vals=(x_star, y_star, s_star), dense=dense)
    end
end


function _cp_to_file(file, params; opt_vals=(), dense=true)
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
        ~(nvars in [5, 8]) && throw(ArgumentError("Invalid file"))
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


function derivatives_to_file(file, deriv, adjoint)
    packed_deriv = (deriv[:dA], deriv[:db], deriv[:dc], deriv[:dx], deriv[:dy], deriv[:ds])
    packed_adjoint = (adjoint[:dA], adjoint[:db], adjoint[:dc], adjoint[:dx], adjoint[:dy], adjoint[:ds])
    writedlm(file, (packed_deriv..., packed_adjoint...))
end


function derivatives_from_file(file)
    file_vec = readlines(file)

    function from_file(file_vec, offset)
        ret = Dict()
        ret[:db] = parse.(Float64, split(file_vec[2+offset], '\t'))
        ret[:dc] = parse.(Float64, split(file_vec[3+offset], '\t'))
        m, n = length(ret[:db]), length(ret[:dc])
        ret[:dA] = reshape(parse.(Float64, split(file_vec[1+offset], '\t')), (m,n))

        ret[:dx] = parse.(Float64, split(file_vec[4+offset], '\t'))
        ret[:dy] = parse.(Float64, split(file_vec[5+offset], '\t'))
        ret[:ds] = parse.(Float64, split(file_vec[6+offset], '\t'))
        return ret
    end

    # derivative, adjoint
    return from_file(file_vec, 0), from_file(file_vec, 6)
end



function _vectorized_index_row_major(i, j)
    if i <= j
        k = div((j-1)*j, 2) + i
    else
        k = div((i-1)*i, 2) + j
    end
    return k
end


function index_map_to_row_major(dim)
    return sortperm(index_map_to_col_major(dim))
end


function index_map_to_col_major(dim)
    index_map = []
    for i=1:dim
        for j=i:dim
            k = _vectorized_index_row_major(i, j)
            push!(index_map, k)
        end
    end
    index_map = sortperm(index_map)
    return index_map
end

function rewrite_sdps_to_row_major!(A, b, cones)
    i = 1
    for cone in cones
        if typeof(cone) <: MOI.PositiveSemidefiniteConeTriangle
            d = MOI.dimension(cone)
            s = cone.side_dimension
            index_map = index_map_to_row_major(cone.side_dimension)
            A[i:(i+d-1), :] .= A[i:(i+d-1), :][index_map, :]
            b[i:(i+d-1)] .= b[i:(i+d-1)][index_map]
        end
        i += MOI.dimension(cone)
    end
end


function rewrite_sdps_to_row_major(A, b, cones)
    i = 1
    A_ = copy(A)
    b_ = copy(b)
    for cone in cones
        if typeof(cone) <: MOI.PositiveSemidefiniteConeTriangle
            d = MOI.dimension(cone)
            s = cone.side_dimension
            index_map = index_map_to_row_major(cone.side_dimension)
            A_[i:(i+d-1), :] .= A[i:(i+d-1), :][index_map, :]
            b_[i:(i+d-1)] .= b[i:(i+d-1)][index_map]
        end
        i += MOI.dimension(cone)
    end
    return A_, b_
end


function rewrite_sdps_to_col_major!(A, b, cones)
    i = 1
    for cone in cones
        if typeof(cone) <: MOI.PositiveSemidefiniteConeTriangle
            d = MOI.dimension(cone)
            s = cone.side_dimension
            index_map = index_map_to_col_major(cone.side_dimension)
            A[i:(i+d-1), :] .= A[i:(i+d-1), :][index_map, :]
            b[i:(i+d-1)] .= b[i:(i+d-1)][index_map]
        end
        i += MOI.dimension(cone)
    end
end


function rewrite_sdps_to_col_major(A, b, cones)
    i = 1
    A_ = copy(A)
    b_ = copy(b)
    for cone in cones
        if typeof(cone) <: MOI.PositiveSemidefiniteConeTriangle
            d = MOI.dimension(cone)
            s = cone.side_dimension
            index_map = index_map_to_col_major(cone.side_dimension)
            A_ .= A[i:(i+d-1), :][index_map, :]
            b_ .= b[i:(i+d-1)][index_map]
        end
        i += MOI.dimension(cone)
    end
    return A_, b_
end
