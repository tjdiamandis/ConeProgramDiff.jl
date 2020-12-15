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

const EXP_CONE_THRESH = 1e-6


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


function MOSD.projection_gradient_on_set(::MOSD.DefaultDistance, v::AbstractVector{T}, ::MOI.ExponentialCone) where {T}
    if v in
end


function _proj_exp_cone(v; dual=false)
    # v has dimension 3
    # Moreau’s decomposition theorem:
    #   v = P_K(v) + P_Ko(v) where K is primal cone and Ko is polar cone = -K*
    #   Dual cone proj: P_K*(v) = -P_Ko(-v)
    v = dual ? -v : v

    if in_exp_cone(v)
        return dual ? zeros(3) : v
    elseif in_exp_cone_dual(v)
        return dual ? -v : zeros(3)
    elseif v[1] <= 0 && v[2] <= 0 #TODO: threshold here??
        return dual ? -[max(-v[1],0); v[2]; 0.0] : [v[1]; 0.0; max(v[3],0)]
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
    # Bisection method for root finding of h(rho) from
    # https://docs.mosek.com/slides/2018/ismp2018/ismp-friberg.pdf
    # h(x) is smooth, strictly increasing, and changes sign
    @views r, s, t = v[1], v[2], v[3]

    # Note: these won't both be infinity (by case 3)
    # TODO: figure out a better replacement for UB
    lb = r > 0 ? 1 - s/r : -1e6
    ub = s > 0 ? r/s : 1e6

    # h(x) = ((x-1)*r + s)/(x^2 - x + 1) * exp(x) - (r - x*s)/(x^2 - x + 1)*exp(-x) - t
    h(x) = (((x-1)*r + s) * exp(x) - (r - x*s)*exp(-x))/(x^2 - x + 1) - t

    rho = find_zero(h, (lb, ub), verbose=false)
    vp = ((rho - 1)*r + s)/(rho^2 - rho + 1) * [rho; 1.0; exp(rho)]
    vd = (r - rho*s)/(rho^2 - rho + 1) * [1.0; 1.0-rho; -exp(-rho)]

    return dual ? -vd : vp
end


# function get_exp_proj_case4(v)
#     # v = [r; s; t]
#     # minimize f = 1/2||x - v||^2
#     # st       g = se^{r/s} - t = 0, s > 0
#     # Solve with primal dual Newton method
#     #   (https://web.stanford.edu/~boyd/papers/pdf/prox_algs.pdf 6.3)
#     ∇g(r,s,t) = [exp(r/s), exp(r/s)*(1-r/s); -1.0]
#     hess_g(r,s,t) = exp(r/s) * [1.0/s  -r/s^2  0;
#                                 -r/s^2 r^2/s^3 0;
#                                 0      0       0]
#     # TODO: init x and λ
#     dxdλ = [I + λ*hess_g(x...) ∇(g...); ∇g(x...)' 0] \ [x - v + λ∇g(x...); x[2]*exp(x[1]/x[2]) - x[3]]
#
#
# end


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
    error("Not Implemented")
    if cone <: MOI.ExponentialCone
        # TODO: Exponential Cone:
        # Grad: http://proceedings.mlr.press/v70/ali17a/ali17a.pdf lemma 3.6
        # Supplement: https://arxiv.org/pdf/1705.00772.pdf

    elseif cone <: MOI.PowerCone
        # TODO: Power cone
        # x = (x, y, z) ∈ R × R × R^n
        # https://link.springer.com/article/10.1007/s00186-015-0514-0
        # or http://www.optimization-online.org/DB_FILE/2014/08/4502.pdf (Theorem 3.1)

    end
end


function MOSD.projection_gradient_on_set(::MOSD.DefaultDistance, v::AbstractVector{T}, ::MOI.Reals) where {T}
    return ones(T, (length(v), length(v)))
end

function MOSD.projection_gradient_on_set(::MOSD.DefaultDistance, v::AbstractVector{T}, ::MOI.Zeros) where {T}
    return zeros(T, (length(v), length(v)))
end
