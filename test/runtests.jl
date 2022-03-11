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
    solver(nlpstop, verbose = true)
    final_nlp_at_x, optimal = nlpstop.current_state, nlpstop.meta.optimal
    @test optimal
    reset!(nlp)
    stats = solver(nlp, verbose = false)
    @test stats.status == :first_order
    reset!(nlp)
end
