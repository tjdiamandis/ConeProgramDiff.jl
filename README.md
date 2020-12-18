# ConeProgramDiff

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://tjdiamandis.github.io/ConeProgramDiff.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://tjdiamandis.github.io/ConeProgramDiff.jl/dev)
[![Build Status](https://github.com/tjdiamandis/ConeProgramDiff.jl/workflows/CI/badge.svg)](https://github.com/tjdiamandis/ConeProgramDiff.jl/actions)
[![Build Status](https://travis-ci.com/tjdiamandis/ConeProgramDiff.jl.svg?branch=master)](https://travis-ci.com/tjdiamandis/ConeProgramDiff.jl)
[![Coverage](https://codecov.io/gh/tjdiamandis/ConeProgramDiff.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/tjdiamandis/ConeProgramDiff.jl)

`ConeProgramDiff.jl` allows the user to compute the derivatives and adjoints of convex cone programs. The package is based on [`diffcp`](https://github.com/cvxgrp/diffcp) and [Differentiating through a cone program](http://web.stanford.edu/~boyd/papers/diff_cone_prog.html).

#### Final project for [18.337](https://mitmath.github.io/18337/)

#### Note: This package is still in active development.

## Cone programs
`ConeProgramDiff.jl` differentiates through a primal-dual cone program pair. The primal problem must be expressed as

```
minimize        c'x
subject to      Ax + s = b
                s in K
```
where  `x` and `s` are variables, `A`, `b` and `c` are the user-supplied problem data, and `K` is a user-defined convex cone. The corresponding dual problem is

```
maximize        -b'y
subject to      A'y + c == 0
                y in K^*
```

with dual variable `y`.


## Usage

We export the function

```julia
solve_and_diff(A, b, c, cone_prod; kwargs)}
```
which takes in the problem parameters `A`, `b`, `c`, and a cone `cone_prod`, which is a vector of the following sets:
* `MOI.Zeros(d)` for the `d`-dimensional zero cone,
* `MOI.Nonnegatives(d)` for the `d`-dimensional positive orthant,
* `MOI.SecondOrderCone(d)` for the `d`-dimensional SOC cone,
* `MOI.PositiveSemidefiniteConeTriangle(d)` for the `d`-dimensional PSD cone (vector version),
* `MOI.ExponentialCone()` for the three-dimensional exponential cones,
* `MOI.DualExponentialCone()` for the three-dimensional dual exponential cones,
* `MOI.PowerCone()` for the three-dimensional power cone, and
* `MOI.DualPowerCone()` for the three-dimensional power cone.
Note that the order does not matter.

Keyword arguments include
* `solver`, which specifies the solver to use for the "forward pass," in which the optimization problem is solved. This defaults to [SCS](https://github.com/cvxgrp/scs). Any JuMP-supported solver can be added by changing `SUPPORTED_SOLVERS` in `cone_solve.jl`
* `use_lsqr`, which is a boolean that indicates if the LSQR method or a dense factorization and backsolve should be used to solve the linear system in the derivative and the adjoint.
* `verbose`, which is a boolean controlling solver output.
* `eps`, which specifies stopping tolerance if supported by the chosen solver.
* `max_iters`, which provides specifies the max iterations for the SCS solver.
* `scale`, which specifies the regularization parameter for SCS. Defaults to 1.0.
* `warm_start`, which provides an optional warm start to the SCS solver.

The method returns a tuple

```julia
x_star, y_star, s_star, derivative, adjoint,
```

The variables  `x_star, y_star, s_star` are the optimal primal, dual, and slack respectively. The function
```julia
derivative(dA, db, dc)
```
takes in a perturbation to the parameters and computes the resulting `dx_star, dy_star, ds_star`.
The function
```julia
adjoint(dx)
```
takes in a perturbation to the solution `x_star` and computes the resulting `dA, db, dc`.

### Additional functionality
In addition to `solve_and_diff`, we provide several utility functions related to cone programming. We export
```julia
project_onto_cone(x, cone_prod)
```
and
```julia
d_project_onto_cone(x, cone_prod)
```
to compute the projection of `x` onto the cone `cone_prod` and its derivative for an arbitrary product cone. The projection is a `Vector` of length `length(x)`, and its derivative is a `Matrix` of size `(length(x), length(x))`.

Additionally, we export utility functions to generate cone programs. Given `dims = (m, n)`, the function
```julia
random_cone_program(dims, cone_prod)
```
generates a random feasible cone program. This function returns a dictionary that includes the parameters `A, b, c` and the solution `x_star, y_star, s_star`. The function
```julia
l1_minimization_program(dims)
```
returns the parameters `A, b, c, cone_prod` of an $\ell_1$ minimization problem, where `size(A) = dims`. The function checks that `dims[1] >= dims[2]` and `rank(A) == dims[2]` so that this problem has a unique solution.

## Future Plans

### Immediate Future
We do not plan to actively maintain this package; instead, we plan to integrate its functionality into the Julia optimization ecosystem. Specifically, we will
* put cone program differentation into [DiffOpt.jl](https://github.com/jump-dev/DiffOpt.jl).
* put projections for the exponential and power cone (and their duals) into [MathOptSetDistances.jl](https://github.com/matbesancon/MathOptSetDistances.jl).

### Longer Term
Ideally, we would like to create something akin to [`cvxpylayers`](https://github.com/cvxgrp/cvxpylayers) for the Julia ecosystem.


## Other
* We compiled mathematical formulations for cone projections and their derivatives in [this document](https://github.com/tjdiamandis/ConeProgramDiff.jl/blob/main/cone_ref.pdf). All mistakes are our own.
* Other code for our final project is in [this repo](https://github.com/csquires/ConeProgramDiff-benchmarking).
