const SUPPORTED_MOSD_INPUT_SETS = Union{
    MOI.Zeros,
    MOI.Nonnegatives,
    MOI.SecondOrderCone,
    MOI.PositiveSemidefiniteConeTriangle,
}

const SUPPORTED_NON_MOSD_INPUT_SETS = Union{
    MOI.ExponentialCone,
    MOI.DualExponentialCone,
    MOI.PowerCone,
    MOI.DualPowerCone
}

const SUPPORTED_MOSD_SETS = Union{
    SUPPORTED_MOSD_INPUT_SETS,
    MOI.Reals
}
const SUPPORTED_NON_MOSD_SETS = SUPPORTED_NON_MOSD_INPUT_SETS

const SUPPORTED_PROJ_SETS = Union{
    SUPPORTED_MOSD_SETS,
    SUPPORTED_NON_MOSD_SETS
}

# Same as cones supported by SCS
const SUPPORTED_INPUT_SETS = Union{
    SUPPORTED_MOSD_INPUT_SETS,
    SUPPORTED_NON_MOSD_INPUT_SETS
}

const EXP_CONE_THRESH = 1e-10
const POW_CONE_THRESH = 1e-14

# Fixes error in MOSD
"""
   unvec_symm(x, dim)
Returns a dim-by-dim symmetric matrix corresponding to `x`.
`x` is a vector of length dim*(dim + 1)/2, corresponding to a symmetric matrix
X = [ X11 X12 ... X1k
       X21 X22 ... X2k
       ...
       Xk1 Xk2 ... Xkk ],
where
vec(X) = (X11, X21, ..., Xk1, X22, X32, ..., Xkk)
"""
function MOSD.unvec_symm(x, dim)
   X = zeros(dim, dim)
   idx = 1
   for i in 1:dim
       for j in 1:i
           # @inbounds X[j,i] = X[i,j] = x[(i-1)*dim-div((i-1)*i, 2)+j]
           X[j,i] = X[i,j] = x[idx]
           idx += 1
       end
   end
   X /= sqrt(2)
   X[LinearAlgebra.diagind(X)] *= sqrt(2)
   return X
end


"""
   vec_symm(X)
Returns a vectorized representation of a symmetric matrix `X`.
`vec(X) = (X11, X21, ..., Xk1, X22, X32, ..., Xkk)`
"""
function MOSD.vec_symm(X)
    X = copy(X)
    X *= sqrt(2)
    X[LinearAlgebra.diagind(X)] .= X[LinearAlgebra.diagind(X)] ./ sqrt(2)
   return X[LinearAlgebra.tril(trues(size(X)))']
end


# Vector <-> matrix
vec_symm = MOSD.vec_symm
unvec_symm = MOSD.unvec_symm


function project_onto_cone(x, cone_prod)
    # TODO: fix to be undefined.
    pi_x = zeros(length(x))
    # pi_x = Vector{typeof(x)}(undef, length(x))
    curr = 1
    for cone in cone_prod
        x_curr = @view(x[curr:curr+MOI.dimension(cone)-1])
        pi_x[curr:curr+MOI.dimension(cone)-1] .= _proj(x_curr, cone)
        curr += MOI.dimension(cone)
    end
    return pi_x
end


function _proj(x, cone)
    if typeof(cone) <: SUPPORTED_MOSD_SETS
        return MOSD.projection_on_set(MOSD.DefaultDistance(), x, cone)
    elseif typeof(cone) <: SUPPORTED_NON_MOSD_SETS
        return _proj_custom(x, cone)
    else
        throw(ArgumentError("Projections for $(typeof(cone)) are not supported"))
    end
end


function _proj_custom(x, cone)
    if typeof(cone) <: MOI.ExponentialCone
        _proj_exp_cone(x; dual=false)
    elseif typeof(cone) <: MOI.DualExponentialCone
        _proj_exp_cone(x; dual=true)

    # TODO: Power cone
    elseif typeof(cone) <: MOI.PowerCone
        error("Not Implemented")
    elseif typeof(cone) <: MOI.DualPowerCone
        error("Not Implemented")
    end
end


function _proj_exp_cone(v; dual=false)
    # v has dimension 3
    # Moreau’s decomposition theorem:
    #   v = P_K(v) + P_Ko(v) where K is primal cone and Ko is polar cone = -K*
    #   Dual cone proj: P_K*(v) = -P_Ko(-v)
    # See section 6.3 of https://web.stanford.edu/~boyd/papers/pdf/prox_algs.pdf
    v = dual ? -v : v

    if in_exp_cone(v)
        return dual ? zeros(3) : v
    elseif in_exp_cone_dual(v)
        return dual ? -v : zeros(3)
    elseif v[1] <= 0 && v[2] <= 0 #TODO: threshold here??
        return dual ? [0.0; -v[2]; max(-v[3],0)] : [v[1]; 0.0; max(v[3],0)]
    else
        return get_exp_proj_case4(v; dual=dual)
    end
end


function in_exp_cone(v)
    return (
        (v[1] <= 0 && abs(v[2]) <= EXP_CONE_THRESH && v[3] >= 0) ||
        (v[2] > 0 && v[2] * exp(v[1] / v[2]) - v[3] <= EXP_CONE_THRESH)
    )
end


function in_exp_cone_dual(v)
    return (
        (abs(v[1]) <= EXP_CONE_THRESH && v[2] <= 0 && v[3] <= 0) ||
        (v[1] > 0 && v[1]*exp(v[2]/v[1]) + ℯ*v[3] <= EXP_CONE_THRESH)
    )
end


function get_exp_proj_case4(v; dual=false)
    # Bisection method for root finding of h(x) from
    # https://docs.mosek.com/slides/2018/ismp2018/ismp-friberg.pdf
    # Thm: h(x) is smooth, strictly increasing, and changes sign
    @views r, s, t = v[1], v[2], v[3]

    # Note: these won't both be infinity (by case 3)
    # TODO: figure out a better replacement for Inf
    lb = r > 0 ? 1 - s/r : -1e6
    ub = s > 0 ? r/s : 1e6

    # h(x) = ((x-1)*r + s)/(x^2 - x + 1) * exp(x) - (r - x*s)/(x^2 - x + 1)*exp(-x) - t
    h(x) = (((x-1)*r + s) * exp(x) - (r - x*s)*exp(-x))/(x^2 - x + 1) - t

    # TODO: anything more efficient than bisection search?
    x = find_zero(h, (lb, ub), verbose=false)
    vp = ((x - 1)*r + s)/(x^2 - x + 1) * [x; 1.0; exp(x)]
    vd = (r - x*s)/(x^2 - x + 1) * [1.0; 1.0-x; -exp(-x)]
# -vd
    return dual ? -v + vp : vp
end


# K = {(x,y,z) | x^a * y^(1-a) >= |z|, x>=0, y>=0}
# K* = {(u,v,w) | (u/a)^a * (v/(1-a))^(1-a) >= |w|, u>=0, v>=0}
function _proj_pow_cone(v, α; dual=false)
    # Prop 2.2 of https://link.springer.com/article/10.1007/s00186-015-0514-0
    # v = dual ? -v : v
    x, y, z = dual ? -v : v

    if in_pow_cone(x, y, z, α)
        # if in power cone
        ret = [x; y; z]
    elseif in_pow_cone_polar(x, y, z, α)
        # if in polar cone Ko = -K*
        ret = zeros(3)
    elseif abs(z) <= 1e-6#POW_CONE_THRESH
        # if not in K, Ko and z = 0
        ret = [max(x,0); max(y,0); 0.0]
    else
        ret = _proj_pow_cone_case4(x, y, z, α)
    end

    return dual ? v + ret : ret
end


function in_pow_cone(x, y, z, α)
    return x >= 0 && y >= 0 && POW_CONE_THRESH + x^α * y^(1-α) >= abs(z)
end


function in_pow_cone_polar(x, y, z, α)
    # -v = -[x;y;z] in K*(α)
    return (x <= 0 && y <=0 &&
        POW_CONE_THRESH + (-x)^α * (-y)^(1-α) >= α^α * (1-α)^(1-α) * abs(z))
end


function _proj_pow_cone_case4(x, y, z, α)
    Phi(x,y,z,α,r) = 0.5*(Phi_prod(x,α,z,r)^α * Phi_prod(y,1-α,z,r)^(1-α)) - r
    Phi_prod(xi,αi,z,r) = (xi + sqrt(xi^2 + 4*αi*r*(abs(z) - r)))

    lb_ub = (0.0+POW_CONE_THRESH, abs(z)-POW_CONE_THRESH)
    # println(Phi(x,y,z,α,lb_ub[1]), Phi(x,y,z,α,lb_ub[2]))
    r = find_zero(r -> Phi(x,y,z,α,r), lb_ub, verbose=false)

    return [0.5*Phi_prod(x,α,z,r); 0.5*Phi_prod(y,1-α,z,r); sign(z)*r]
end


# TODO: can probably write this code more efficiently? lots of repeats
function d_project_onto_cone(x, cone_prod)
    # TODO: fix to be undefined.
    ret = Vector{SparseMatrixCSC}(undef, length(cone_prod))
    # pi_x = Vector{typeof(x)}(undef, length(x))
    curr = 1
    for (ind, cone) in enumerate(cone_prod)
        ret[ind] = sparse(_d_proj(@view(x[curr:curr+MOI.dimension(cone)-1]), cone))
        curr += MOI.dimension(cone)
    end
    return blockdiag(ret...)
end


function _d_proj(x, cone)
    if typeof(cone) <: SUPPORTED_MOSD_SETS
        return MOSD.projection_gradient_on_set(MOSD.DefaultDistance(), x, cone)
    elseif typeof(cone) <: SUPPORTED_NON_MOSD_SETS
        return _d_proj_custom(x, cone)
    else
        throw(ArgumentError("Cone type $(typeof(cone)) is not supported."))
    end
end


function _d_proj_custom(x, cone)
    if typeof(cone) <: MOI.ExponentialCone
        _d_proj_exp_cone(x; dual=false)
    elseif typeof(cone) <: MOI.DualExponentialCone
        _d_proj_exp_cone(x; dual=true)

    # TODO: Power cone
    elseif typeof(cone) <: MOI.PowerCone
        error("Not Implemented")
    elseif typeof(cone) <: MOI.DualPowerCone
        error("Not Implemented")
    end
end


function _d_proj_exp_cone(v; dual=false)
    # Lemma 3.6 https://arxiv.org/pdf/1705.00772.pdf
    # J_PK*(v) = I - J_PK(-v)
    v = dual ? -v : v
    Ip(z) = z >= 0 ? 1 : 0

    if in_exp_cone(v)
        return dual ? zeros(3,3) : Matrix{Float64}(I, 3, 3)
    elseif in_exp_cone_dual(v)
        return dual ? Matrix{Float64}(I, 3, 3) : zeros(3,3)
    elseif v[1] <= 0 && v[2] <= 0 #TODO: threshold here??
        vp = diagm([1; Ip(v[2]); Ip(v[3])])
        return dual ? I - vp : vp
    else
        # Have to look at KKT conditions of projection problem
        return get_exp_d_proj_case4(v; dual=dual)
    end
end


function get_exp_d_proj_case4(v; dual=false)
    Pv = get_exp_proj_case4(v, dual=false)
    nu = Pv[3] - v[3]
    z1,z2,z3 = Pv
    abs(Pv[2]) <= 1e-8 && error("Numerical error in projection onto exp cone")
    rs = z1/z2
    exp_rs = exp(rs)

    mat = inv([
        1+nu*exp_rs/z2     -nu*exp_rs*rs/z2       0     exp_rs;
        -nu*exp_rs*rs/z2   1+nu*exp_rs*rs^2/z2    0     (1-rs)*exp_rs;
        0                  0                      1     -1
        exp_rs             (1-rs)*exp_rs          -1    0
    ])
    return @view(mat[1:3,1:3])
end


function _d_proj_pow_cone(v, α; dual=false)
    # Thm 3.1 of https://link.springer.com/article/10.1007/s00186-015-0514-0
    x, y, z = dual ? -v : v

    if in_pow_cone(x, y, z, α)
        # if in power cone
        ret = Matrix{Float64}(I, 3, 3)
    elseif in_pow_cone_polar(x, y, z, α)
        # if in polar cone Ko = -K*
        ret = zeros(3,3)
    elseif abs(z) <= POW_CONE_THRESH
        # if not in K, Ko and z = 0
        ret = _d_proj_pow_cone_case3(x, y, z, α)
    else
        ret = _d_proj_pow_cone_case4(x, y, z, α)
    end

    return dual ? Matrix{Float64}(I, 3, 3) - ret : ret
end


function _d_proj_pow_cone_case3(x, y, z, α)
    v = [x; y]
    I(t) = t > 0 ? 1 : 0
    αs = [α; 1-α]

    if sum(αs[v .> 0]) > sum(αs[v .< 0])
        d = 1
    elseif sum(αs[v .> 0]) < sum(αs[v .< 0])
        d = 0
    else
        num = reduce(*, (-v[v .< 0]).^αs[v .< 0])
        denom = reduce(*, v[v .> 0]^αs[v .> 0]) * reduce(*, αs[v .< 0]^αs[v .< 0])
        d = 1/((num/denom)^2 + 1)
    end
    return diagm([I(x), I(y), d])
end


function _d_proj_pow_cone_case4(x, y, z, α)
    Phi(x,y,z,α,r) = 0.5*(Phi_prod(x,α,z,r)^α * Phi_prod(y,1-α,z,r)^(1-α)) - r
    Phi_prod(xi,αi,z,r) = (xi + sqrt(xi^2 + 4*αi*r*(abs(z) - r)))

    lb_ub = (0.0+POW_CONE_THRESH, abs(z)-POW_CONE_THRESH)
    r = find_zero(r -> Phi(x,y,z,α,r), lb_ub, verbose=false)

    za = abs(z)
    gx = sqrt(x^2 + 4*α*r*(za - r))
    gy = sqrt(y^2 + 4*(1-α)*r*(za - r))
    fx = 0.5*(x + gx)
    fy = 0.5*(y + gy)

    β = 1-α
    T_neg = (α*x/gx + β*y/gy)
    L = 2*(za - r) / (za + (za - 2r) * T_neg)
    T = -T_neg
    J_ii(w, γ, g) = 0.5 + w/2g + γ^2*(za - 2r)*r*L/g^2
    J_ij = α*β*(za - 2r)*r*L/(gx*gy)
    J = [
        J_ii(x, α, gx)      J_ij                sign(z)*α*r*L/gx;
        J_ij                J_ii(y, β, gy)      sign(z)*β*r*L/gy;
        sign(z)*α*r*L/gx    sign(z)*β*r*L/gy    r/za*(1+T*L)
    ]
    return J
end


function MOSD.projection_gradient_on_set(::MOSD.DefaultDistance, v::AbstractVector{T}, ::MOI.Reals) where {T}
    return Matrix{Float64}(I, (length(v), length(v)))
end

function MOSD.projection_gradient_on_set(::MOSD.DefaultDistance, v::AbstractVector{T}, ::MOI.Zeros) where {T}
    return zeros(T, (length(v), length(v)))
end


function MOSD. projection_gradient_on_set(::MOSD.DefaultDistance, v::AbstractVector{T}, ::MOI.PositiveSemidefiniteConeTriangle) where {T}
    dim = isqrt(2*length(v))
    X = unvec_symm(v, dim)
    λ, U = LinearAlgebra.eigen(X)
    D = LinearAlgebra.Diagonal(max.(λ, 0))

end
