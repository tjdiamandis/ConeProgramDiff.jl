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
