# Vector <-> matrix
vec_symm = MOSD.vec_symm
unvec_symm = MOSD.unvec_symm


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
