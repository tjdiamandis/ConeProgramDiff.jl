module ConeProgramDiff

using LinearAlgebra, SparseArrays                       #stdlib
using JuMP, MathOptInterface, MathOptSetDistances       #Opt utils
using IterativeSolvers, Roots                           #Eq solvers
using SCS, ECOS, Hypatia                                #Opt solvers
using Random

const MOI = MathOptInterface
const MOSD = MathOptSetDistances

include("utils.jl")
include("cones.jl")
include("cone_solve.jl")
include("generate_program.jl")

export solve_and_diff
export project_onto_cone, d_project_onto_cone
export random_cone_program, l1_minimization_program

end


# TODO: Check out https://github.com/Jutho/LinearMaps.jl for projections
