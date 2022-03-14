using Pkg; Pkg.activate("")
using CUTEst, Dates, NLPModels, SolverBenchmark

using ARCTR, JSOSolvers, NLPModelsIpopt

nmax = 10
problems = readlines("list_problems_$nmax.dat")
cutest_problems = (CUTEstModel(p) for p in problems)

max_time = 120.0 #20 minutes
max_ev = typemax(Int)
max_iter = typemax(Int)
tol = 1e-5

solvers = Dict(
  :ARCqKOp => nlp -> ARCqKOp(
    nlp,
    verbose = false,
    atol = tol,
    rtol = tol,
    max_time = max_time,
    max_iter = max_iter,
  ),
  :ST_ARCOp => nlp -> ST_ARCOp(
    nlp,
    verbose = false,
    atol = tol,
    rtol = tol,
    max_time = max_time,
    max_iter = max_iter,
  ),
  :ST_TROp => nlp -> ST_TROp(
    nlp,
    verbose = false,
    atol = tol,
    rtol = tol,
    max_time = max_time,
    max_iter = max_iter,
  ),
)
stats = bmark_solvers(solvers, cutest_problems)

using JLD2
@save "$(today())_$(prod(String.(keys(solvers)) .* "_"))cutest_$(string(length(problems)))_$(nmax).jld2" stats
