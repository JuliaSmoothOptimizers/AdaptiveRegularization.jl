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
    println(nbsolver, "  ", solver)
    eval(solver)(nlpstop, verbose = true)
    final_nlp_at_x, optimal = nlpstop.current_state, nlpstop.meta.optimal
    @test optimal
    reset!(nlp)
    stats = eval(solver)(nlp, verbose = false)
    @test stats.status == :first_order
    reset!(nlp)
end
