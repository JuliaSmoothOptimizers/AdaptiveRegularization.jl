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
# nlp = MathProgNLPModel(woods(), name = "Woods");
# nlpatx = NLPAtX(nlp.meta.x0)
# nlpstop = NLPStopping(nlp, Stopping.unconstrained, nlpatx)
#
# final_nlp_at_x, optimal = TRARC2(nlp, nlpstop, verbose = true)
#
# printstyled("on vérifie avec la version PAS stopping \n", color = :cyan)
# nlp2 = MathProgNLPModel(woods(), name = "Woods");
# x,f, norm_grad, norm_grad_opt, niter, calls, OK, status = TRARC(nlp2, nlp.meta.x0, TrustRegion(10.0), Combi(hessian_dense, PDataLDLt, solve_modelTRDiag, preprocessLDLt, decreaseFact, Tparam()), verbose = true)

# global nbsolver = 0
# for solver in ALL_solvers
#     global nbsolver += 1
#     println(nbsolver,"  ",solver)
#     (x, f, gNorm, iter, optimal, tired, status) = solver(nlp, verbose = true)
#     @printf("%-15s  %8d  %9.2e  %7.1e  %5d  %5d  %6d  %s\n",
#             nlp.meta.name, nlp.meta.nvar, f, gNorm,
#             nlp.counters.neval_obj, nlp.counters.neval_grad,
#             nlp.counters.neval_hprod, status)    #stats = run_solver(solver, model, verbose=false)
#     #@test (all([stats...] .>= 0))
#     @test optimal
#     reset!(nlp)
# end

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
    #@test (all([stats...] .>= 0))
    @test optimal
    reset!(nlp)
end
