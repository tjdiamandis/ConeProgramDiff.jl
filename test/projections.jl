include("utils.jl")

function test_proj_zero()
    Random.seed!(0)
    n = 100
    for _ in 1:10
        x = randn(n)
        @assert zeros(n) ≈ ConeProgramDiff._proj(x, MOI.Zeros(n))
        @assert x ≈ ConeProgramDiff._proj(x, MOI.dual_set(MOI.Zeros(n)))
    end
    return true
end

function test_proj_pos()
    Random.seed!(0)
    n = 100
    for _ in 1:10
        x = randn(n)
        @assert max.(x,zeros(n)) ≈ ConeProgramDiff._proj(x, MOI.Nonnegatives(n))
        @assert max.(x,zeros(n)) ≈ ConeProgramDiff._proj(x, MOI.dual_set(MOI.Nonnegatives(n)))
    end
    return true
end


function test_proj_soc()
    Random.seed!(0)
    n = 100
    for _ in 1:10
        x = randn(n)
        model = Model()
        set_optimizer(model, optimizer_with_attributes(SCS.Optimizer, "eps" => 1e-8, "max_iters" => 10000, "verbose" => 0))
        @variable(model, z[1:n])
        @variable(model, t)
        @objective(model, Min, t)
        @constraint(model, sum((x-z).^2) <= t)
        @constraint(model, z in MOI.SecondOrderCone(n))
        optimize!(model)

        p = ConeProgramDiff._proj(x, MOI.SecondOrderCone(n))
        @assert p ≈ value.(z)
        @assert p ≈ ConeProgramDiff._proj(x, MOI.dual_set(MOI.SecondOrderCone(n)))
    end
    return true
end


# TODO: Fix PSD cone -- Triangle vs square matrix
function test_proj_psd()
    Random.seed!(0)
    n = 10
    for _ in 1
        x = randn(n,n)
        x = x + x'
        x_vec = x[LinearAlgebra.tril(trues(size(x)))']
        model = Model()
        set_optimizer(model, optimizer_with_attributes(SCS.Optimizer, "eps" => 1e-8, "max_iters" => 10000, "verbose" => 0))
        @variable(model, z[1:n,1:n])
        @variable(model, t)
        @objective(model, Min, t)
        @constraint(model, sum((x-z).^2) <= t)
        @constraint(model, z in PSDCone())
        optimize!(model)

        z_star = value.(z)[LinearAlgebra.tril(trues(size(x)))']
        p = ConeProgramDiff._proj(x_vec, MOI.PositiveSemidefiniteConeTriangle(n))
        if !isapprox(p, z_star)
            error("Error in PSD Cone Projection:\np=$p\nz=$z_star")
        end
        @assert p ≈ ConeProgramDiff._proj(x_vec, MOI.dual_set(MOI.PositiveSemidefiniteConeTriangle(n)))
    end
    return true
end

# Test Exponential Cone
function det_case_exp_cone(v; dual=false)
    v = dual ? -v : v
    if ConeProgramDiff.in_exp_cone(v)
        return 1
    elseif ConeProgramDiff.in_exp_cone_dual(v)
        return 2
    elseif v[1] <= 0 && v[2] <= 0 #TODO: threshold here??
        return 3
    else
        return 4
    end
end

function test_proj_exp(tol)
    function _test_proj_exp_cone_help(x, tol; dual=false)
        cone = dual ? MOI.DualExponentialCone() : MOI.ExponentialCone()
        model = Model()
        set_optimizer(model, optimizer_with_attributes(SCS.Optimizer, "eps" => 1e-8, "max_iters" => 10000, "verbose" => 0))
        @variable(model, z[1:3])
        @variable(model, t)
        @objective(model, Min, t)
        @constraint(model, sum((x-z).^2) <= t)
        @constraint(model, z in cone)
        optimize!(model)
        z_star = value.(z)
        @assert isapprox(ConeProgramDiff._proj_exp_cone(x,  dual=dual), z_star, atol=tol)
    end

    Random.seed!(0)
    n = 3
    case_p = zeros(4)
    case_d = zeros(4)
    for _ in 1:100
        # x = ConeProgramDiff._proj_exp_cone(x) + randn(3)
        x = randn(3)
        # println(x)
        case_p[det_case_exp_cone(x; dual=false)] += 1
        _test_proj_exp_cone_help(x, tol; dual=false)

        case_d[det_case_exp_cone(x; dual=true)] += 1
        _test_proj_exp_cone_help(x, tol; dual=true)
    end
    @assert all(case_p .> 0) && all(case_d .> 0)
    return true
end


function det_case_pow_cone(v, α; dual=false)
    v = dual ? -v : v
    x, y, z = v

    if ConeProgramDiff.in_pow_cone(x, y, z, α)
        # if in power cone
        return 1
    elseif ConeProgramDiff.in_pow_cone_polar(x, y, z, α)
        # if in polar cone Ko = -K*
        return 2
    elseif abs(z) <= 1e-6
        # if not in K, Ko and z = 0
        return 3
    else
        return 4
    end
end

function test_proj_pow_cone(tol)
    function _test_proj_pow_cone_help(x, α, tol; dual=false)
        cone = dual ? MOI.DualPowerCone(α) : MOI.PowerCone(α)
        model = Model()
        set_optimizer(model, optimizer_with_attributes(
            SCS.Optimizer, "eps" => 1e-8, "max_iters" => 10000,
            "verbose" => 0, "scale" => 10.0))
        @variable(model, z[1:3])
        @variable(model, t)
        @objective(model, Min, t)
        @constraint(model, sum((x-z).^2) <= t)
        @constraint(model, z in cone)
        optimize!(model)
        z_star = value.(z)
        xp = ConeProgramDiff._proj_pow_cone(x, α, dual=dual)
        if !isapprox(xp, z_star, atol=tol)
            error("x = $x,\nz* = $z_star,\nxp = $xp,\ncase = $(det_case_pow_cone(x, α; dual=dual))")
        end
    end
    Random.seed!(0)
    n = 3
    case_p = zeros(4)
    case_d = zeros(4)
    for _ in 1:100
        x = randn(3)
        α = rand()*0.6 + 0.2
        # α = rand()

        # Need to get some into case 3
        if rand(1:10) == 1
            x[3] = 0
        end
        case_p[det_case_pow_cone(x, α; dual=false)] += 1
        _test_proj_pow_cone_help(x, α, tol; dual=false)

        case_d[det_case_pow_cone(x, α; dual=true)] += 1
        _test_proj_pow_cone_help(x, α, tol; dual=true)
    end

    @assert all(case_p .> 0) && all(case_d .> 0)
    return true
end


function test_project_onto_cone()
    function _test_proj_all_cones(x, cone_prod; dual=false)
        cones = dual ? [MOI.dual_set(c) for c in cone_prod] : cone_prod
        xp = project_onto_cone(x, cones)
        offset = 0
        for cone in cones
            inds = 1+offset:offset+MOI.dimension(cone)
            @assert xp[inds] ≈ ConeProgramDiff._proj(x[inds], cone)
            offset += MOI.dimension(cone)
        end
    end

    for _ in 1:10
        cone_prod = get_random_product_cone(10)
        N = sum([MOI.dimension(c) for c in cone_prod])
        local x = randn(N)
        _test_proj_all_cones(x, cone_prod; dual=false)
        _test_proj_all_cones(x, cone_prod; dual=true)
    end
    return true
end



# Test sets
# ------------------------------------------------------------------------
@testset "Projections Test (Zero, Pos, SOC, PSD)" begin
    @test test_proj_zero()
    @test test_proj_pos()
    @test test_proj_soc()
    @test test_proj_psd()
end

@testset "Projection Test: Exp Cone" begin
    @test test_proj_exp(1e-7)
end

@testset "Projection Test: Pow Cone" begin
    @test test_proj_pow_cone(1e-6)
end

@testset "Projection: project_onto_cone()" begin
    test_project_onto_cone()
end
