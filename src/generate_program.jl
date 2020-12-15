using DelimitedFiles

# Generates a non-trivial random cone program
# min c'x
#  st Ax + s = b
#          s in K
# K is a convex cone (product of given cones)
function random_cone_program(dims, cone_prod)
    m,n = dims

    # Break into ortho parts in K and in K*
    z = randn(m)
    s_star = project_onto_cone(z, cone_prod)
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


function l1_minimization_program(dims; solve=false)
    m,n = dims
    m >= n || throw(ArgumentError("m < n"))
    A = randn(dims)
    b = randn(m)
    @assert rank(A) == n

    c_ = [zeros(n); ones(m)]
    A_ = -[
        -A                          Matrix{Float64}(I, m, m);
        A                           Matrix{Float64}(I, m, m);
        Matrix{Float64}(I, n, n)    zeros(n, m);
        ones(n)'                    zeros(m)'
    ]
    b_ = [b; -b; zeros(n); -1.0]
    # -Ax + b == s ∈ K
    cones = [MOI.Nonnegatives(2m+n), MOI.Zeros(1)]

    if solve
        model = Model()
        set_optimizer(model, optimizer_with_attributes(SCS.Optimizer,
            "eps" => 1e-10, "verbose" => 0))
        @variable(model, x[1:10])
        @variable(model, t[1:20])
        @objective(model, Min, sum(t))
        @constraint(model, 0 .<= A*x - b + t)
        @constraint(model, 0 .<= t - A*x + b)
        @constraint(model, x .>= 0)
        @constraint(model, sum(x) == 1.0)
        optimize!(model)
        return A_, b_, c_, cones, model
    else
        return A_, b_, c_, cones
    end

end

#
# A, b, c, cone_prod, model = l1_minimization_program((20, 10))
# m,n = 20,10
# xs = [value(x) for x in all_variables(model)][1:n]
# ts = value.(all_variables(model)[n+1:end])
# opt_val = sum(ts)
# sum(abs.(A[1:m,1:n]*xs - b[1:m]))
# all(A*[xs;ts]  - b .>= -1e10)
# (A*[xs;ts]  - b)[end] <= 1e10
# dot(c,[xs;ts])
#
# solve_opt_problem(A, b, c, cone_prod, nothing, SCS.Optimizer)
# # solve_opt_problem(A, b, c, [MOI.Nonnegatives(size(A)[1])], nothing, SCS.Optimizer)
# m,n = size(A)
# using Convex
# x = Variable(n)
# s = Variable(m)
# p = minimize(dot(x,c), A*x + s == b, s >= 0)
# solve!(p, () -> SCS.Optimizer())
#
# model = Model()
# optimizer_constructor = optimizer_with_attributes(SCS.Optimizer, "eps" => 1e-10, "verbose" => 1)
# set_optimizer(model, optimizer_constructor)
# @variable(model, x[1:size(A)[2]])
# @variable(model, s[1:size(A)[1]])
# @objective(model, Min, dot(c,x))
# @constraint(model, b - A*x .>= s)
# @constraint(model, s[end:end] in MOI.Zeros(1))
# @constraint(model, s[1:end-1] .>= 0)
# optimize!(model)
# optimizer_constructor
# MOI.get(owner_model(x[1]), MOI.VariablePrimal(1), x[1])
# typeof(backend(model)) <: SCS.Optimizer
# MOI.optimize!(backend(model))
#
# solution = SCS_solve(m, n, A, b, c, f=0, l=);
