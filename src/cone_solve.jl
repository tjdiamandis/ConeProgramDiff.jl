#
# diffcp.ZERO for the zero cone,
# diffcp.POS for the positive orthant,
# diffcp.SOC for a product of SOC cones,
# diffcp.PSD for a product of PSD cones, and
# diffcp.EXP for a product of exponential cones.

# TODO: Maybe say zero_cone instead of zero
# TODO: cone_dict(ConeProgramDiff.zero_cone => dim)
@enum Cone zero pos


function Base.in(vec::Vector, cone::Cone)
    if cone == zero
        return all(vec .== 0)
    end
end

function add_constraint!(model, s, cone)
    if cone == zero
        @constraint(model, s .== 0)
    else

    end
end


model = Model(Hypatia.Optimizer)
@variable(model, x[1:10])
@constaint(A*x + s .== c)

# TODO: Type function
function solve_and_diff(A, b, c, cone_dict, warm_start=nothing, solver=nothing) T <: Float
    m,n = size(A)
    # Cone dict check: valid cones and dimension match
    @assert all([c in CONES for (c, n) in cone_dict])
    @assert sum([n for (c, n) in cone_dict]) == m
    # dimension check
    @assert length(b) == m
    @assert length(c) == n

    # Solver select
    # TODO: Potentially allow more cones for Hypatia
    if solver == "Hypatia"
        model = Model(Hypatia.Optimizer)
    elseif solver == "ECOS"
        model = Model(ECOS.Optimizer)
    elseif solver == "SCS" || isnothing(solver)
        model = Model(SCS.Optimizer)
    else
        throw(ArgumentError("Invalid solver"))
    end

    # Model
    @variable(model, x[1:n])
    @variable(model, s[1:m])
    @objective(model, c'*x)
    @constaint(A*x + s .== c)

    curr = 1
    for (cone, dim) in cone_dict
        #TODO: overload in for cones?
        add_constraint(model, s[], cone)
        @constraint(model, s[curr:curr+dim] in cone)
        curr += dim
    end

    optimize!(model)
    xstar = value.(x)

end


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
@code_warntype @constraint(model, vcat(hypo, [Q[i, j] for i in 1:2 for j in 1:i]...) in MOI.RootDetConeTriangle(2))
@variable(model, hypo)

# solve
optimize!(model)
termination_status(model)
objective_value(model)
value.(x)

([2,1] in MOI.RootDetConeTriangle(2))
