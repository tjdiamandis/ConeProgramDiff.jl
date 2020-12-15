module ConeProgramDiff

using LinearAlgebra, SparseArrays                       #stdlib
using JuMP, MathOptInterface, MathOptSetDistances       #Opt utils
using IterativeSolvers
using SCS, ECOS, Hypatia                                #Solvers

const MOI = MathOptInterface
const MOSD = MathOptSetDistances

include("utils.jl")
include("cones.jl")
include("cone_solve.jl")
include("generate_program.jl")

export random_cone_program
export solve_and_diff
export project_onto_cone, d_project_onto_cone

end
