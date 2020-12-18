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
    x_star, y_star, s_star, sol = solve_opt_problem(A, b, c, cone_prod, warm_start, optimizer)

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
    return x_star, y_star, s_star, pushforward, pullback, sol
end

cone_order = Dict(
    MOI.Zeros=>1,
    MOI.Nonnegatives=>2,
    MOI.SecondOrderCone=>3,
    MOI.PositiveSemidefiniteConeTriangle=>4,
    MOI.ExponentialCone=>5,
    MOI.DualExponentialCone=>6,
    MOI.PowerCone=>7
)

function reorder_opt_problem_scs(A, b, cones)
    # sort the cones so that they're in the correct order
    sort_ixs = sortperm(cones, by=cone->cone_order[typeof(cone)])
    sorted_cones = cones[sort_ixs]

    f = 0  # the number of zero cones
    l = 0  # the number of linear cones
    q = Int[]  # the array of SOC sizes
    s = Int[]  # the array of SDC sizes
    ep = 0  # the number of exponential cones
    ed = 0  # the number of dual exponential cones
    p = Float64[]  # the array of power cone parameters

    # create the parameters for SCS
    for cone in sorted_cones
        if typeof(cone) <: MOI.Zeros
            f += MOI.dimension(cone)
        elseif typeof(cone) <: MOI.Nonnegatives
            l += MOI.dimension(cone)
        elseif typeof(cone) <: MOI.SecondOrderCone
            push!(q, MOI.dimension(cone))
        elseif typeof(cone) <: MOI.PositiveSemidefiniteConeTriangle
            push!(s, MOI.dimension(cone))
        elseif typeof(cone) <: MOI.ExponentialCone
            ep += 1
        elseif typeof(cone) <: MOI.DualExponentialCone
            ep += 1
        elseif typeof(cone) <: MOI.PowerCone
            push!(p, MOI.dimension(a))
        elseif typeof(cone) <: MOI.DualPowerCone
            push!(p, -MOI.dimension(a))
        end
    end
    cone_dict = Dict(:f=>f, :l=>l, :q=>q, :s=>s, :ep=>ep, :ed=>ed, :p=>p)

    original_ixs = []
    start_ix = 1
    for cone in cones
        end_ix = start_ix + MOI.dimension(cone) - 1
        push!(original_ixs, start_ix:end_ix)
        start_ix = end_ix + 1
    end

    index_map = [original_ixs[i] for i in sort_ixs]
    index_map = collect(Iterators.flatten(index_map))

    return A[index_map, :], b[index_map], cone_dict
end


function solve_opt_problem(A, b, c, cone_prod, warm_start, optimizer_factory)
    # TODO: potentially play with the scale depending on the cones?
    #       e.g. see https://github.com/jump-dev/Convex.jl/issues/104
    #       Essentially want to have the primal and dual res down at the same rate
    #           primal slower -> inc. scale
    #           dual slower -> dec. scale
    m,n = size(A)

    A_reordered, b_reordered, cone_dict = reorder_opt_problem_scs(A, b, cone_prod)
    f, l, q, s = cone_dict[:f], cone_dict[:l], cone_dict[:q], cone_dict[:s]
    ep, ed, p = cone_dict[:ep], cone_dict[:ed], cone_dict[:p]
    sol = SCS_solve(SCS.DirectSolver, m, n, A_reordered, b_reordered, c, f, l, q, s, ep, ed, p)
    println(sol)

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
    return value.(x), y, value.(s), sol
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
