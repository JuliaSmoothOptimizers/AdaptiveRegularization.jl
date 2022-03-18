using Pkg; Pkg.activate("")
using CUTEst, Dates, NLPModels, SolverBenchmark

using ARCTR, JSOSolvers, NLPModelsIpopt

nmax = 10
problems = readlines("list_problems_$nmax.dat")
cutest_problems = (CUTEstModel(p) for p in problems)

max_time = 3600.0
max_ev = typemax(Int)
max_iter = typemax(Int)
atol = 1e-5
rtol = 1e-6

solvers = Dict(
  :ARCqKOp => nlp -> ARCqKOp(
    nlp,
    verbose = false,
    atol = atol,
    rtol = rtol,
    max_time = max_time,
    max_iter = max_iter,
  ),
  :ST_TROp => nlp -> ST_TROp(
    nlp,
    verbose = false,
    atol = atol,
    rtol = rtol,
    max_time = max_time,
    max_iter = max_iter,
  ),
)
stats = bmark_solvers(solvers, cutest_problems)

using JLD2
@save "$(today())_$(prod(String.(keys(solvers)) .* "_"))cutest_$(string(length(problems)))_$(nmax).jld2" stats
