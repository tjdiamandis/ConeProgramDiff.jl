using ConeProgramDiff


using MathOptInterface
using SCS
using Random
using LinearAlgebra

const MOI = MathOptInterface


function check_kkt(A, b, c, cone_prod, x, y, s; tol=1e-8)
    try
        @assert isapprox(A*x + s, b,atol=tol)
        @assert isapprox(A'*y + c, zeros(length(c)),atol=tol)
        @assert isapprox(s, project_onto_cone(s, cone_prod),atol=tol)
        @assert isapprox(y, project_onto_cone(y, [MOI.dual_set(c) for c in cone_prod]),atol=tol)
        @assert isapprox(s'*y, 0, atol=tol)
    catch e
        error("KKT conditions failed:
               \nprimal constraint: $(norm(A*x + s - b))
               \ndual constraint: $(norm(A'*y +c))
               \nK error: $(norm(s - project_onto_cone(s, cone_prod)))
               \nK* error: $(norm(y - project_onto_cone(y, [MOI.dual_set(c) for c in cone_prod])))
               \nCS cond: $(norm(s'*y))"
        )
    end
end


function check_adjoint(tol)
    Random.seed!(0)
    dims = (15, 10)
    A, b, c, cone_prod = l1_minimization_program(dims; solve=false)
    x, y, s, pf, pb = solve_and_diff(A, b, c, cone_prod; eps=1e-10, max_iters=10000)
    check_kkt(A, b, c, cone_prod, x, y, s, tol=1e-8)


    # Check adjoint
    pstar = dot(x, c)
    dA, db, dc = pb(c)
    del = 1e-6
    x_pert, _, _, _, _ = solve_and_diff(A+del*dA, b+del*db, c+del*dc, cone_prod;
                                                  eps=1e-10, max_iters=10000)
    pstar_pert = dot(x_pert, c)
    df = pstar_pert - pstar
    dp = del * (sum(dA .* dA) + db'*db + dc'*dc)
    if !isapprox(df, dp, atol=tol)
        error("Adjoint mismatch: $(abs(df-dp))\n\tdf = $df\n\tdp = $dp")
    end
    return true
end


function check_derivative()
    Random.seed!(0)
    dims = (15, 10)
    A, b, c, cone_prod = l1_minimization_program(dims; solve=false)
    x, y, s, pf, pb = solve_and_diff(A, b, c, cone_prod; eps=1e-10, max_iters=10000)
    check_kkt(A, b, c, cone_prod, x, y, s, tol=1e-6)

    # Check derivative
    del = 1e-6
    dA, db, dc = del*randn(size(A)), del*randn(size(b)), del*randn(size(c))
    dx, dy, ds = pf(dA, db, dc)
    x_pert, y_pert, s_pert, _, _ = solve_and_diff(A+dA, b+db, c+dc, cone_prod;
                                                  eps=1e-10, max_iters=10000)

    # TODO: tols are high here -- should investigate
    @assert isapprox(x_pert - x, dx, atol=1e-4)
    @assert isapprox(y_pert - y, dy, atol=1e-4)
    @assert isapprox(s_pert - s, ds, atol=2e-4)
    return true
end


@testset "solve_and_diff Tests" begin
    @test check_adjoint(2e-5)
    @test check_derivative()
end
