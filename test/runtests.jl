using ARCTR
using LinearAlgebra
using Test
using Printf
using NLPModels
using NLPModelsJuMP
using JuMP

# first test the LDLt function
include("../src/Utilities/testLDLt.jl")

# test all solvers with the well known Woods test function
include("rosenbrock.jl")
nlp = MathProgNLPModel(extrosnb(2), name = "extrosnb");

global nbsolver = 0
for solver in ALL_solvers
    global nbsolver += 1
    println(nbsolver,"  ",solver)
    (x, f, gNorm, iter, optimal, tired, status) = solver(nlp, verbose = true)
    @printf("%-15s  %8d  %9.2e  %7.1e  %5d  %5d  %6d  %s\n",
            nlp.meta.name, nlp.meta.nvar, f, gNorm,
            nlp.counters.neval_obj, nlp.counters.neval_grad,
            nlp.counters.neval_hprod, status)    #stats = run_solver(solver, model, verbose=false)
    #@test (all([stats...] .>= 0))
    @test optimal
    reset!(nlp)
end
