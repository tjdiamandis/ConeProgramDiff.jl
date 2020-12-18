using ConeProgramDiff
using Test
using LinearAlgebra, SparseArrays                       #stdlib
using JuMP, MathOptInterface, MathOptSetDistances       #Opt utils
using IterativeSolvers, Roots                           #Eq solvers
using SCS, ECOS, Hypatia                                #Opt solvers
using Random

const MOI = MathOptInterface
const MOSD = MathOptSetDistances

# @testset "ConeProgramDiff.jl" begin
#     # Write your tests here.
# end

@testset "Projections" begin
    include("projections.jl")
end

@testset "D_Projections" begin
    include("d_projections.jl")
end

@testset "Solve and diff" begin
    include("solve_and_diff.jl")
end
