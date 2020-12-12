const SUPPORTED_SOLVERS = Dict(
    "ECOS" => ECOS.Optimizer,
    "Hypatia" => Hypatia.Optimizer,
    "SCS" => SCS.Optimizer
)

# TODO: Type function params A,b,c
function solve_and_diff(
    A, b, c, cone_prod::Vector{T}; warm_start=nothing, solver="SCS"
) where {T <: MOI.AbstractVectorSet}
    m,n = size(A)

    # Arguments check
    (length(b) == m && length(c) == n) || throw(DimensionMismatch("Mismatch between A, b, c"))
    all([typeof(c) <: SUPPORTED_INPUT_SETS for c in cone_prod]) || throw(ArgumentError("Unsupported cones"))
    all([MOI.dimension(c) > 0 for c in cone_prod]) || throw(ArgumentError("Cones must have nonzero dimension"))
    sum([MOI.dimension(c) for c in cone_prod]) == m || throw(DimensionMismatch("Mismatch between dimension of c and cone"))
    if solver in keys(SUPPORTED_SOLVERS)
        optimizer = SUPPORTED_SOLVERS[solver]
        # TODO: check cone_prod for each solver
    else
        throw(ArgumentError("Invalid solver"))
    end

    return _solve_and_diff(A, b, c, cone_prod, warm_start, optimizer)
end

function _solve_and_diff(A, b, c, cone_prod, warm_start, optimizer)
    m,n = size(A)
    x_star, y_star, s_star = solve_opt_problem(A, b, c, cone_prod, warm_start, optimizer)

    Q = spzeros(m+n+1,m+n+1)
    Q[1:n,n+1:n+m]      = A'
    Q[n+1:n+m,1:n]      = -A
    Q[n+1:n+m,m+n+1]    = b
    Q[m+n+1,n+1:n+m]    = -b'
    Q[1:n,m+n+1]        = c
    Q[m+n+1,1:n]        = -c'

    DπKdual_v = d_project_onto_cone(y_star - s_star, cone_prod)
    println(DπKdual_v)

    # TODO: extend to dy and ds
    function pullback(dx, dy, ds)
        u, v, w = (x_star, y_star - s_star, 1.0) #def z
        dz = [
            dx;
            DπKdual_v'*(dy + ds) - ds;
            -x_star' * dx - y_star' * dy - s_star' * ds
        ]
        Pi_z = pi_z(u, v, w, cone_prod)

        M = (Q - I) * dpi(u, v, w, cone_prod) + I
        g = lsqr(-M', dz)
        dQ = g * Pi_z'

        dA = dQ[1:n,n+1:n+m]'  - dQ[n+1:n+m,1:n]    # dQ12' - dQ21
        db = dQ[n+1:n+m,m+n+1] - dQ[m+n+1,n+1:n+m]  # dQ23 - dQ32'
        dc = dQ[1:n,m+n+1]     - dQ[m+n+1,1:n]      # dQ13 - dQ31'
        return dA, db, dc
    end

    function pushforward(dA, db, dc)
        u, v, w = (x_star, y_star - s_star, 1.0)
        dQ = spzeros(n+m+1, n+m+1)
        dQ[1:n,n+1:n+m]      = dA'
        dQ[n+1:n+m,1:n]      = -dA
        dQ[n+1:n+m,m+n+1]    = db
        dQ[m+n+1,n+1:n+m]    = -db'
        dQ[1:n,m+n+1]        = dc
        dQ[m+n+1,1:n]        = -dc'
        M = (Q - I) * dpi(u, v, w, cone_prod) + I
        g = dQ*pi_z(u, v, w, cone_prod)
        dz = lsqr(-M, g)
        du, dv, dw = dz[1:n], dz[n+1:n+m], dz[n+m+1]
        dx = du .- dw*x_star
        dy = DπKdual_v*dv - dw*y_star
        ds = DπKdual_v*dv - dv - dw*s_star
        return dx, dy, ds
    end
    return x_star, y_star, s_star, pushforward, pullback
end

function solve_opt_problem(A, b, c, cone_prod, warm_start, optimizer_factory)
    m,n = size(A)
    model = Model(optimizer_factory)
    @variable(model, x[1:n])
    @variable(model, s[1:m])
    @objective(model, Min, c'*x)
    con = @constraint(model, A*x + s .== b)
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
    return value.(x), dual.(con), value.(s)
end

function pi_z(u, v, w, cone_prod)
    return vcat(
        u,
        project_onto_cone(v, [MOI.dual_set(c) for c in cone_prod]),
        max(w, 0.0)
    )
end

function dpi(u, v, w, cone_prod)
    return blockdiag(
        sparse(Matrix{Float64}(I, length(u), length(u))),
        d_project_onto_cone(v, [MOI.dual_set(c) for c in cone_prod]),
        sparse(Matrix{Float64}(I, 1,1))
    )
end
