using ARCTR
using LinearAlgebra
using Test
using Printf
using NLPModels
using Stopping

using ADNLPModels
using OptimizationProblems.ADNLPProblems

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
