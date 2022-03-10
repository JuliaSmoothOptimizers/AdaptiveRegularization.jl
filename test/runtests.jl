using ARCTR
using LinearAlgebra
using Test
using Printf
using NLPModels
using Stopping

using ADNLPModels
using OptimizationProblems.ADNLPProblems

include("testLDLt.jl")
# include("testAbs.jl") #  MethodError: no method matching ARCSpectral(::ADNLPModel{Float64, Vector{Float64}}; verbose=true)
# include("testMA57_0.jl")
# include("testMA97.jl")

global nbsolver = 0
for solver in ALL_solvers
    global nbsolver += 1
    nlp = extrosnb(n = 2)
    nlpstop = NLPStopping(nlp)
    println(nbsolver,"  ",solver)
    final_nlp_at_x, optimal = solver(nlp, nlpstop, verbose = true)
    @test optimal
    reset!(nlp)
end
