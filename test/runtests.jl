# This package
using ARCTR
# stdlib
using LinearAlgebra, SparseArrays, Test
# JSO
using ADNLPModels, NLPModels, OptimizationProblems.ADNLPProblems, SolverTest
# Stopping
using Stopping

@testset "Testing NLP solvers" begin
    @testset "$name" for name in keys(ARCTR.solvers_const)
        solver = eval(name)
        unconstrained_nlp(solver)
        multiprecision_nlp(solver, :unc)
    end
end

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

if VERSION >= v"1.7.0"
    include("allocation_test.jl")
    include("allocation_test_main.jl")
end
