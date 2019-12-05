using ARCTR
using LinearAlgebra
using Test
using Printf
using NLPModels
using Stopping


# test all solvers with the well known Rosenbrock test function
include("rosenbrock.jl")

global nbsolver = 0
for solver in ALL_solvers
    global nbsolver += 1
    nlp = ADNLPModel(extrosnb, [-1.0, -1.0], name = "rosenbrock")
    nlpatx = NLPAtX([-1.0, -1.0])
    nlpstop = NLPStopping(nlp, Stopping.unconstrained, nlpatx)
    println(nbsolver,"  ",solver)
    final_nlp_at_x, optimal = solver(nlp, nlpstop, verbose = true)
    @test optimal
    reset!(nlp)
end
