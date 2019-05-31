using ARCTR
using LinearAlgebra
using Test
using Printf
using NLPModels
using NLPModelsJuMP
using JuMP
using State
using Stopping


# first test the LDLt function
# include("../src/Utilities/testLDLt.jl")

# test all solvers with the well known Woods test function
include("woods.jl")

global nbsolver = 0
for solver in ALL_solvers_stopping
    global nbsolver += 1
    nlp = MathProgNLPModel(woods(), name = "Woods");
    nlpatx = NLPAtX(nlp.meta.x0)
    nlpstop = NLPStopping(nlp, Stopping.unconstrained, nlpatx)
    println(nbsolver,"  ",solver)
    final_nlp_at_x, optimal = solver(nlp, nlpstop, verbose = true)
    @printf("%-15s  %8d  %9.2e  %7.1e  %5d  %5d  %6d  %s\n",
            nlp.meta.name, nlp.meta.nvar, final_nlp_at_x.fx, norm(final_nlp_at_x.gx),
            nlp.counters.neval_obj, nlp.counters.neval_grad,
            nlp.counters.neval_hprod, optimal)    #stats = run_solver(solver, model, verbose=false)
    @test optimal
    reset!(nlp)
end
