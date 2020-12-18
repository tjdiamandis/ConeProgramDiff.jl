const SUPPORTED_SOLVERS = Dict(
    "ECOS" => ECOS.Optimizer,
    "Hypatia" => Hypatia.Optimizer,
    "SCS" => SCS.Optimizer
)

# TODO: Type function params A,b,c
function solve_and_diff(
    A::AbstractMatrix{S}, b::AbstractVector{S}, c::Vector{S}, cone_prod::Vector{T}; kwargs...
) where {T <: MOI.AbstractVectorSet, S <: AbstractFloat}
    kwargs = Dict(kwargs)
    m,n = size(A)

    # Arguments check
    (length(b) == m && length(c) == n) || throw(DimensionMismatch("Mismatch between A, b, c"))
    all([typeof(c) <: SUPPORTED_INPUT_SETS for c in cone_prod]) || throw(ArgumentError("Unsupported cones"))
    all([MOI.dimension(c) > 0 for c in cone_prod]) || throw(ArgumentError("Cones must have nonzero dimension"))
    sum([MOI.dimension(c) for c in cone_prod]) == m || throw(DimensionMismatch("Mismatch between dimension of c and cone"))
    if haskey(kwargs, :solver) && kwargs[:solver] in keys(SUPPORTED_SOLVERS)
        optimizer = SUPPORTED_SOLVERS[solver]
        # TODO: check cone_prod for each solver
    elseif haskey(kwargs, :solver)
        throw(ArgumentError("Invalid solver"))
    else
        optimizer = SCS.Optimizer
    end

    warm_start = haskey(kwargs, :warm_start) ? kwargs[:warm_start] : nothing
    if haskey(kwargs, :method) && kwars[:method] == "dense"
        use_lsqr = false
    else
        use_lsqr = true
    end

    return _solve_and_diff(A, b, c, cone_prod, warm_start, optimizer, use_lsqr)
end


function _solve_and_diff(A, b, c, cone_prod, warm_start, optimizer, use_lsqr)
    m,n = size(A)
    typeof(A) <: SparseMatrixCSC && dropzeros!(A)
    x_star, y_star, s_star = solve_opt_problem(A, b, c, cone_prod, warm_start, optimizer)

    Q = spzeros(m+n+1,m+n+1)
    Q[1:n,n+1:n+m]      .= A'
    Q[n+1:n+m,1:n]      .= -A
    Q[n+1:n+m,m+n+1]    .= b
    Q[m+n+1,n+1:n+m]    .= -b
    Q[1:n,m+n+1]        .= c
    Q[m+n+1,1:n]        .= -c

    u, v, w = (x_star, y_star - s_star, 1.0)
    M = (Q - I) * dpi_z(u, v, w, cone_prod) + I

    function pullback(dx)
        # Assume that dy = ds = 0
        dz = [
            dx;
            zeros(length(y_star));
            -x_star' * dx
        ]
        Pi_z = pi_z(u, v, w, cone_prod)

        g = use_lsqr ? lsqr(-M', dz, atol=1e-10, btol=1e-10) : -M' \ dz

        # TODO: make more efficient by only pulling needed cols/rows
        dQ = g * pi_z(u, v, w, cone_prod)'
        dA = dQ[1:n,n+1:n+m]'  - dQ[n+1:n+m,1:n]    # dQ12' - dQ21
        db = dQ[n+1:n+m,m+n+1] - dQ[m+n+1,n+1:n+m]  # dQ23 - dQ32'
        dc = dQ[1:n,m+n+1]     - dQ[m+n+1,1:n]      # dQ13 - dQ31'
        return dA, db, dc
    end

    function pushforward(dA, db, dc)
        # TODO: Could we compute a matrix vector product instead of full matrix?
        DπKdual_v = d_project_onto_cone(v, [MOI.dual_set(c) for c in cone_prod])

        dQ = spzeros(n+m+1, n+m+1)
        dQ[1:n,n+1:n+m]      .= dA'
        dQ[n+1:n+m,1:n]      .= -dA
        dQ[n+1:n+m,m+n+1]    .= db
        dQ[m+n+1,n+1:n+m]    .= -db #dont use tranpose bc of broadcast
        dQ[1:n,m+n+1]        .= dc
        dQ[m+n+1,1:n]        .= -dc

        g = dQ*pi_z(u, v, w, cone_prod)

        dz = use_lsqr ? lsqr(-M, g, atol=1e-10, btol=1e-10) : -M \ g
        @views du, dv, dw = dz[1:n], dz[n+1:n+m], dz[n+m+1]
        dx = du .- dw*x_star
        dy = DπKdual_v*dv - dw*y_star
        ds = DπKdual_v*dv - dv - dw*s_star
        return dx, dy, ds
    end
    return x_star, y_star, s_star, pushforward, pullback
end


function solve_opt_problem(A, b, c, cone_prod, warm_start, optimizer_factory)
    # TODO: potentially play with the scale depending on the cones?
    #       e.g. see https://github.com/jump-dev/Convex.jl/issues/104
    #       Essentially want to have the primal and dual res down at the same rate
    #           primal slower -> inc. scale
    #           dual slower -> dec. scale
    m,n = size(A)
    model = Model()
    set_optimizer(model, optimizer_with_attributes(
        SCS.Optimizer, "eps" => 1e-10, "max_iters" => 100000, "verbose" => 0))
    @variable(model, x[1:n])
    @variable(model, s[1:m])
    @objective(model, Min, c'*x)
    con = @constraint(model, A*x + s .== b)
    curr = 1
    for cone in cone_prod
        @constraint(model, s[curr:curr+MOI.dimension(cone)-1] in cone)
        curr += MOI.dimension(cone)
    end
    optimize!(model)
    if !(termination_status(model) in [MOI.OPTIMAL, MOI.ALMOST_OPTIMAL])
        error("Model not solved correctly.
               Termination status: $(termination_status(model))")
    end
    if termination_status(model) != MOI.OPTIMAL
        @warn "Maximum iterations reached. Set 'max_iters' to increase accuracy"
    end

    # primal_status(model)
    # dual_status(model)
    # TODO: for some reason this fails when b is a sparse vector?
    #    Maybe interface directly with SCS?
    y = isapprox(A'* dual.(con) + c, zeros(length(c)), atol=1e-6) ? dual.(con) : -dual.(con)
    return value.(x), y, value.(s)
end


function pi_z(u, v, w, cone_prod; dual=false)
    cone = dual ? [MOI.dual_set(c) for c in cone_prod] : cone_prod
    return vcat(
        u,
        project_onto_cone(v, cone),
        max(w, 0.0)
    )
end


function dpi_z(u, v, w, cone_prod; dual=false)
    cone = dual ? [MOI.dual_set(c) for c in cone_prod] : cone_prod
    return blockdiag(
        sparse(Matrix{Float64}(I, length(u), length(u))),
        d_project_onto_cone(v, cone),
        sparse(Matrix{Float64}(I, 1,1))
    )
end
