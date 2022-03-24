using Pkg;
Pkg.activate("");
using CUTEst, Dates, NLPModels, SolverBenchmark

using ARCTR, JSOSolvers, NLPModelsIpopt

nmax = 1000000
problems = readlines("list_problems_$nmax.dat")
cutest_problems = (CUTEstModel(p) for p in problems)

max_time = 3600.0
max_ev = typemax(Int)
max_iter = typemax(Int)
atol = 1e-5
rtol = 1e-6

solvers = Dict(
    :ARCqKOp10201 =>
        nlp -> ARCqKOp(
            nlp,
            verbose = false,
            atol = atol,
            rtol = rtol,
            max_time = max_time,
            max_iter = max_iter,
            shifts = 10.0 .^ (collect(-10.0:1.0:20.0)),
        ),
    :ARCqKOp10202 =>
        nlp -> ARCqKOp(
            nlp,
            verbose = false,
            atol = atol,
            rtol = rtol,
            max_time = max_time,
            max_iter = max_iter,
            shifts = 10.0 .^ (collect(-10.0:2.0:20.0)),
        ),
    :ARCqKOp10203 =>
        nlp -> ARCqKOp(
            nlp,
            verbose = false,
            atol = atol,
            rtol = rtol,
            max_time = max_time,
            max_iter = max_iter,
            shifts = 10.0 .^ (collect(-10.0:3.0:20.0)),
        ),
    :ARCqKOp10204 =>
        nlp -> ARCqKOp(
            nlp,
            verbose = false,
            atol = atol,
            rtol = rtol,
            max_time = max_time,
            max_iter = max_iter,
            shifts = 10.0 .^ (collect(-10.0:4.0:20.0)),
        ),
)
stats = bmark_solvers(solvers, cutest_problems)
