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

const SUPPORTED_SOLVERS = Dict(
    "ECOS" => ECOS.Optimizer,
    "Hypatia" => Hypatia.Optimizer,
    "SCS" => SCS.Optimizer
)


# TODO: Type function params A,b,c
function solve_and_diff(
    A, b, c, cone_prod::Vector{T}; warm_start=nothing, solver="SCS"
) where {T <: SUPPORTED_INPUT_SETS}
    m,n = size(A)

    # Arguments check
    @assert length(b) == m
    @assert length(c) == n
    @assert all([MOI.dimension(c) > 0 for c in cone_prod])
    @assert sum([MOI.dimension(c) for c in cone_prod]) == m
    if solver in keys(SUPPORTED_SOLVERS)
        optimizer = SUPPORTED_SOLVERS[solver]
        # TODO: check cone_prod for each solver
    else
        throw(ArgumentError("Invalid solver"))
    end

    return _solve_and_diff(A, b, c, cone_prod, warm_start, optimizer)
end

function _solve_and_diff(
    A, b, c, cone_prod::Vector{T}, warm_start, optimizer
) where {T <: SUPPORTED_INPUT_SETS}
    m,n = size(A)
    x_star, y_star, s_star = solve_opt_problem(A, b, c, cone_prod, warm_start, optimizer)

    # backward pass
    # TODO: extend to dy and ds
    function pullback(dx)
        u, v, w = (x_star, y_star - s_star, 1.0) #def z
        dz = [dx, zeros(m), -x' * dx]
        Pi_z = pi_z(u, w, v)
        M = (Q - I) * dpi(u, v, w, cone_prod) + I
        g = lsqr(-M', dz)
        dQ = g' * Pi_z

        dA = dQ[1:n,n+1:n+m]'  + dQ[n+1:n+m,1:n]    # dQ12' - dQ21
        db = dQ[n+1:n+m,m+n+1] + dQ[m+n+1,n+1:n+m]  # dQ23 - dQ32'
        dc = dQ[1:n,m+n+1]     + dQ[m+n+1,1:n]'     # dQ13 - dQ31'
        return dA, db, dc
    end

    function push_forward(dA, db, dc)
        u, v, w = (x_star, y_star - s_star, 1.0)
        dQ = spzeros(n+m+1, n+m+1)
        # Set dA, db, dc
        # M = (Q - I) * dpi(u, v, w, cone_prod) + I
        # g = dQ*pi_z(u, v, w)
        # dz = lsqr(-M, g)
        du, dv, dw = dz[1:n], dz[n+1:n+m], dz[n+m+1]
        dx = du .- dw*x
        dy = nothing        #proj onto dual cone TODO
        ds = nothing        #proj onto dual cone TODO
        return dx, dy, ds
    end
    return x_star, y_star, s_star, pushforward, pullback
end

function solve_opt_problem(A, b, c, cone_prod, warm_start, optimizer_factory)
    model = Model(optimizer_factory)
    @variable(model, x[1:n])
    @variable(model, s[1:m])
    @objective(model, Min, c'*x)
    @constraint(model, A*x + s .== b)
    curr = 1
    for cone in cone_prod
        @constraint(model, s[curr:curr+cone.dimension-1] in cone)
        curr += cone.dimension
    end
    optimize!(model)
    if ~(termination_status(model) == MOI.OPTIMAL)
        error("Model not solved correctly.
               Termination status: $(termination_status(model))")
    end

    # primal_status(model)
    # dual_status(model)
    return value.(x), dual.(model), value.(s)
end

function pi_z(u, v, w)
    return vcat(
        u,
        project_onto_cone(v, [MOI.dual_set(c) for c in cone_prod]),
        max(w, 0.0)
    )
end

function dpi(u, v, w, cone_prod)
    return blockdiag(
        Matrix{Float64}(I, length(u), length(u)),
        _d_proj(c, cone_prod),
        Matrix{Float64}(I, 1,1)
    )
end

function project_onto_cone(x, cone_prod)
    # TODO: fix to be undefined.
    pi_x = zeros(length(x))
    # pi_x = Vector{typeof(x)}(undef, length(x))
    curr = 1
    for cone in cone_prod
        x_curr = @view(x[curr:curr+MOI.dimension(cone)-1])
        pi_x[curr:curr+MOI.dimension(cone)-1] .= _proj(x_curr, cone)
        curr += cone.dimension
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


function _d_proj(x, cone)
    if typeof(cone) <: SUPPORTED_MOSD_SETS
        return MOSD.projection_gradient_on_set(MOSD.DefaultDistance(), x, cone)
    elseif typeof(cone) <: SUPPORTED_NON_MOSD_SETS
        return _d_proj_custom(x, cone)
    else
        throw(ArgumentError("Cone is not of supported type."))
    end
end


function _proj_custom(x, cone::SUPPORTED_NON_MOSD_SETS)
    error("Not Implemented")
    # TODO: Exponential Cone:
    # Projection: https://web.stanford.edu/~boyd/papers/pdf/prox_algs.pdf 6.3
    # TODO: Power cone
    # https://link.springer.com/article/10.1007/s00186-015-0514-0
end


function _d_proj_custom(x, cone::SUPPORTED_NON_MOSD_CONES)
    error("Not Implemented")
    # TODO: Exponential Cone:
    # Grad: http://proceedings.mlr.press/v70/ali17a/ali17a.pdf lemma 3.6
    # TODO: Power cone
    # https://link.springer.com/article/10.1007/s00186-015-0514-0
end
