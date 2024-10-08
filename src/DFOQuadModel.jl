module DFOQuadModel

using LagrangePoly
using JuMP
using SumOfSquares
using CSDP
using DynamicPolynomials
using SemialgebraicSets

include("main.jl")
include("subproblem.jl")

export DFTR
export optimize!
export solve_subproblem

end
