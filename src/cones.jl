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
    error("Not Implemented")
    # TODO: Exponential Cone:
    # Projection: https://web.stanford.edu/~boyd/papers/pdf/prox_algs.pdf 6.3
    # TODO: Power cone
    # https://link.springer.com/article/10.1007/s00186-015-0514-0
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
    error("Not Implemented")
    # TODO: Exponential Cone:
    # Grad: http://proceedings.mlr.press/v70/ali17a/ali17a.pdf lemma 3.6
    # TODO: Power cone
    # https://link.springer.com/article/10.1007/s00186-015-0514-0
end

function MOSD.projection_gradient_on_set(::MOSD.DefaultDistance, v::AbstractVector{T}, ::MOI.Reals) where {T}
    return ones(T, (length(v), length(v)))
end

function MOSD.projection_gradient_on_set(::MOSD.DefaultDistance, v::AbstractVector{T}, ::MOI.Zeros) where {T}
    return zeros(T, (length(v), length(v)))
end
