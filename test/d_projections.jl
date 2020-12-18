include("utils.jl")

function _test_d_proj(cone, tol)
    n = MOI.dimension(cone)
    x = randn(n)
    dx = 1e-6 * randn(n)
    proj_x = ConeProgramDiff._proj(x, cone)
    proj_xdx = ConeProgramDiff._proj(x+dx, cone)
    dproj_finite_diff = proj_xdx - proj_x

    dproj_test = ConeProgramDiff._d_proj(x, cone) * dx
    if !isapprox(dproj_finite_diff, dproj_test, atol=tol)
        println("Error!")
        error("Error for x = $x, dx = $dx, $cone")
    end
end


function test_d_proj(cone::MOI.AbstractVectorSet; tol=1e-8)
    Random.seed!(0)
    for _ in 1:100
        _test_d_proj(cone, tol)
        _test_d_proj(MOI.dual_set(cone), tol)
    end
end


function test_d_proj(n::Int; tol=1e-8)
    cones = [
        MOI.Zeros(n),
        MOI.Nonnegatives(n),
        MOI.SecondOrderCone(n),
        MOI.PositiveSemidefiniteConeTriangle(n)
    ]
    for cone in cones
        # println("Testing $cone")
        test_d_proj(cone; tol=tol)
    end
    return true
end


function test_d_proj_exp(tol)
    function _helper(x, dx, tol; dual)
        proj_x = ConeProgramDiff._proj_exp_cone(x,  dual=dual)
        proj_xdx = ConeProgramDiff._proj_exp_cone(x+dx,  dual=dual)
        dproj_finite_diff = proj_xdx - proj_x
        dproj_test = ConeProgramDiff._d_proj_exp_cone(x; dual=dual) * dx
        # println(dproj_finite_diff)
        # println(dproj_test)
        @assert isapprox(dproj_finite_diff, dproj_test, atol=tol)
    end

    Random.seed!(0)
    case_p = zeros(4)
    case_d = zeros(4)
    for _ in 1:100
        x = randn(3)
        dx = 1e-6 * randn(3)
        case_p[det_case_exp_cone(x; dual=false)] += 1
        # println(det_case_exp_cone(x; dual=false))
        _helper(x, dx, tol; dual=false)

        case_d[det_case_exp_cone(x; dual=true)] += 1
        # println(det_case_exp_cone(x; dual=true))
        # println(x)
        _helper(x, dx, tol; dual=true)
    end
    @assert all(case_p .> 0) && all(case_d .> 0)
    return true
end


function test_d_proj_pow(tol)
    function _helper(x, α, dx, tol; dual)
        # println("$x \t $dx \t $(x+dx)")
        proj_x = ConeProgramDiff._proj_pow_cone(x, α; dual=dual)
        proj_xdx = ConeProgramDiff._proj_pow_cone(x+dx, α; dual=dual)
        dproj_finite_diff = proj_xdx - proj_x
        dproj_test = ConeProgramDiff._d_proj_pow_cone(x, α; dual=dual) * dx
        # println(dproj_finite_diff)
        # println(dproj_test)
        @assert isapprox(dproj_finite_diff, dproj_test, atol=tol)
    end

    Random.seed!(0)
    case_p = zeros(4)
    case_d = zeros(4)
    for _ in 1:100
        x = randn(3)
        dx = 1e-6 * randn(3)
        α = rand()*0.6 + 0.2

        # Need to get some into case 3
        if rand(1:10) == 1
            x[3] = 0
        end

        case_p[det_case_pow_cone(x, α; dual=false)] += 1
        _helper(x, α, dx, tol; dual=false)

        case_d[det_case_pow_cone(x, α; dual=true)] += 1
        _helper(x, α, dx, tol; dual=true)
    end
    @assert all(case_p .> 0) && all(case_d .> 0)
    return true
end


function test_d_project_onto_cone()
    function _test_d_proj_all_cones(x, cone_prod; tol=1e-8, dual=false)
        cones = dual ? [MOI.dual_set(c) for c in cone_prod] : cone_prod
        dx = 1e-7 * randn(length(x))
        Dpi = d_project_onto_cone(x, cones)
        proj_x = ConeProgramDiff.project_onto_cone(x, cones)
        proj_xdx = ConeProgramDiff.project_onto_cone(x+dx, cones)
        dproj_finite_diff = proj_xdx - proj_x

        dproj_test = Dpi * dx
        @assert isapprox(dproj_finite_diff, dproj_test, atol=tol)
        # if ~isapprox(dproj_finite_diff, dproj_test, atol=tol)
        #     println("Test failed:")
        #     println("determinant: $(det(Dpi))")
        #     inds = abs.(dproj_finite_diff - dproj_test) .>= tol
        #     println((dproj_finite_diff - dproj_test)[inds])
        #     println()
        # end
    end

    for _ in 1:10
        cone_prod = get_random_product_cone(10)
        N = sum([MOI.dimension(c) for c in cone_prod])
        local x = randn(N)
        _test_d_proj_all_cones(x, cone_prod; dual=false)
        _test_d_proj_all_cones(x, cone_prod; dual=true)
    end
    return true
end


@testset "D_Projections Test (Zero, Pos, SOC, PSD)" begin
    @test test_d_proj(10)
end

@testset "D_Projection Test: Exp Cone" begin
    @test test_d_proj_exp(1e-5)
end

@testset "D_Projection Test: Pow Cone" begin
    @test test_d_proj_pow(1e-6)
end

@testset "D_Projection: d_project_onto_cone()" begin
    test_d_project_onto_cone()
end
