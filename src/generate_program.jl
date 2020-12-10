using DelimitedFiles

# Generates a non-trivial random cone program
# min c'x
#  st Ax + s = b
#          s in K
# K is a convex cone (product of given cones)
function random_cone_program(dims, cone_dict)
    m,n = dims

    # Break into ortho parts in K and in K*
    z = randn(m)
    s_star = _pi(z, cone_dict)
    y_star = z - s_star
    A = sparse(randn(dims))
    x_star = randn(n)
    b = A*x_star + s_star
    c = -A' * y_star

    # TODO: return as a dict instead of tuple?
    return Dict(:A => A, :b => b, :c => c, :x_star => x_star,
                :y_star => y_star, :s_star => s_star)
    # return (A, b, c), (x_star, y_star, s_star)
end


function cp_to_file(file, params, vals=())
    writedlm(file, (params..., vals...))
end

function cp_from_file(file)
    file_vec = readlines("test.txt")
    nvars = length(file_vec)
    if  ~(nvars in [3, 6])
        throw(ArgumentError("Invalid file"))
    end

    A_flat = parse.(Float64, split(file_vec[1], '\t'))
    b = parse.(Float64, split(file_vec[2], '\t'))
    c = parse.(Float64, split(file_vec[3], '\t'))

    m, n = length(b), length(c)
    A = reshape(A_flat, (m,n))

    if nvars == 6
        x_star = parse.(Float64, split(file_vec[4], '\t'))
        y_star = parse.(Float64, split(file_vec[5], '\t'))
        s_star = parse.(Float64, split(file_vec[6], '\t'))
    else
        x_star, y_star, s_star = nothing, nothing, nothing
    end
    return Dict(:A => A, :b => b, :c => c, :x_star => x_star,
                :y_star => y_star, :s_star => s_star)

end
