#
# diffcp.ZERO for the zero cone,
# diffcp.POS for the positive orthant,
# diffcp.SOC for a product of SOC cones,
# diffcp.PSD for a product of PSD cones, and
# diffcp.EXP for a product of exponential cones.

SUPPORTED_MOSD_SETS = Union{MOI.SecondOrderCone, MOI.Nonnegatives}

model = Model(Hypatia.Optimizer)
@variable(model, x[1:10])
@constaint(A*x + s .== c)

# TODO: Type function
function solve_and_diff(A, b, c, cone_dict; warm_start=nothing, solver="SCS")  # T <: Float
    m,n = size(A)
    # Cone dict check: valid cones and dimension match
    # TODO: check that cone is in supported cones
    # @assert all([c in CONES for (c, n) in cone_dict])
    # check dimension matches
    @assert all([c.dimension > 0 for c in cone_dict])
    @assert sum([c.dimension for c in cone_dict]) == m
    # dimension check
    @assert length(b) == m
    @assert length(c) == n

    # Solver select
    # TODO: Potentially allow more cones for Hypatia
    if solver == "Hypatia"
        model = Model(Hypatia.Optimizer)
    elseif solver == "ECOS"
        model = Model(ECOS.Optimizer)
    elseif solver == "SCS"
        model = Model(SCS.Optimizer)
    else
        throw(ArgumentError("Invalid solver"))
    end

    # Model
    @variable(model, x[1:n])
    @variable(model, s[1:m])
    @objective(model, Max, c'*x)
    @constraint(model, A*x + s .== b)

    curr = 1
    for cone in cone_dict
        @constraint(model, s[curr:curr+cone.dimension-1] in cone)
        curr += cone.dimension
    end

    optimize!(model)
    xstar = value.(x)

    # backward pass

    function pullback(dx)
        dz = [dx, zeros(m), -x' * dx]
        M = nothing  # TODO: page 5, requires projection gradient
        g = nothing  # TODO: linear system solve of M and dz
    end
    # dz


    # dQ - use sparse matrices

    # overall projection
    #   - projections onto each cone
end

function _pi(x, cone_dict)
    # TODO: fix to be undefined.
    pi_x = zeros(length(x))
    curr = 1
    for cone in cone_dict
        a = _proj(x[curr:curr+cone.dimension-1], cone)
        pi_x[curr:curr+cone.dimension-1] .= a
        curr += cone.dimension
    end
    return pi_x
end

function _proj(x, cone)
    if typeof(cone) <: SUPPORTED_MOSD_SETS
        return MOSD.projection_on_set(MOSD.DefaultDistance(), x, cone)
    else
        throw(ArgumentError("Cone is not of supported type."))
    end
end

function _proj_grad(x, cone)
    if typeof(cone) <: SUPPORTED_MOSD_SETS
        return MOSD.projection_gradient_on_set(MOSD.DefaultDistance(), x, cone)
    else
        throw(ArgumentError("Cone is not of supported type."))
    end
end

_pi(rand(4), [MOI.Nonnegatives(2), MOI.SecondOrderCone(2)])

function rand_cone_prog(m, n)
    A = rand(m, n)
    b = rand(m)
    c = rand(n)
    return A, b, c
end
A, b, c = rand_cone_prog(3, 4)
solve_and_diff(A, b, c, [MOI.Nonnegatives(3)], solver="Hypatia")


function test_solve()
    items  = [:Gold, :Silver, :Bronze]
    values = Dict(:Gold => 5.0,  :Silver => 3.0,  :Bronze => 1.0)
    weight = Dict(:Gold => 2.0,  :Silver => 1.5,  :Bronze => 0.3)

    model = Model(SCS.Optimizer)
    @variable(model, 0 <= take[items] <= 1)  # Define a variable for each item
    @objective(model, Max, sum(values[item] * take[item] for item in items))
    @constraint(model, sum(weight[item] * take[item] for item in items) <= 3)
    optimize!(model)
    println(value.(take))
    return value
end

items  = [:Gold, :Silver, :Bronze]
vals = Dict(:Gold => 5.0,  :Silver => 3.0,  :Bronze => 1.0)
weight = Dict(:Gold => 2.0,  :Silver => 1.5,  :Bronze => 0.3)
model = Model(SCS.Optimizer)
@variable(model, 0 <= take[items] <= 1)  # Define a variable for each item
@objective(model, Max, sum(vals[item] * take[item] for item in items))
@variable(model, x[1:10])
x
@constraint(model, x in zero)
take
# setup model
V = rand(2, 3)

model = Model(Hypatia.Optimizer)
@variable(model, x[1:3] >= 0)
@constraint(model, sum(x) == 5)
@variable(model, hypo)
@objective(model, Max, hypo)
Q = V * diagm(x) * V'
@constraint(model, vcat(hypo, [Q[i, j] for i in 1:2 for j in 1:i]...) in MOI.RootDetConeTriangle(2))
@variable(model, hypo)

# solve
optimize!(model)
termination_status(model)
objective_value(model)
value.(x)

using SCS
using LinearAlgebra
using Hypatia
using JuMP
using MathOptInterface
using MathOptSetDistances
const MOI = MathOptInterface
const MOSD = MathOptSetDistances

A = rand(3) .- .5
MOSD.projection_on_set(MOSD.DefaultDistance(), A, MOI.Nonnegatives(3))
eye4 = [1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]  # symmetrical PSD triangle format
MOSD.projection_on_set(MOSD.DefaultDistance(), eye4, MOI.PositiveSemidefiniteConeTriangle(4))
A = rand(4)
MOSD.projection_on_set(MOSD.DefaultDistance(), A, MOI.SecondOrderCone(4))
MOSD.projection_gradient_on_set

Dict("a"=>1)
Dict(MOI.Zeros=>3, MOI.Reals=>2, MOI.Zeros=>2)
[(MOI.Zeros(3)), (MOI.Reals(2))]
