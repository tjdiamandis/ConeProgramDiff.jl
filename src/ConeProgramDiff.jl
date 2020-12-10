module ConeProgramDiff

using LinearAlgebra, SparseArrays       #stdlib
using JuMP, MathOptInterface            #Opt utils
using SCS, ECOS, Hypatia                #Solvers

const MOI = MathOptInterface

# include("cone_solve.jl")
include("generate_program.jl")

# export project_onto_cone

end
