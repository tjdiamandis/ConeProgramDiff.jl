using ConeProgramDiff
using Test

using MathOptInterface
const MOI = MathOptInterface
# const MOSD = MathOptSetDistances

@testset "ConeProgramDiff.jl" begin
    # Write your tests here.
end

@testset "Projections" begin
    include("projections.jl")
end
