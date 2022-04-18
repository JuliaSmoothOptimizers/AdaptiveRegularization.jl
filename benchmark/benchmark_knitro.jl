using Pkg;
Pkg.activate("knitro");
using CUTEst, Dates, NLPModels, SolverBenchmark
using NLPModelsKnitro

nmax = 10
problems = readlines("list_problems_$nmax.dat")
cutest_problems = (CUTEstModel(p) for p in problems)

max_time = 3600.0
max_ev = typemax(Int)
max_iter = typemax(Int)
atol = 1e-5
rtol = 1e-6

solvers = Dict(
    :knitro =>
        nlp -> knitro(
            nlp,
            opttol = rtol,
            opttol_abs = atol,
            maxfevals = max_ev,
            maxit = max_iter,
            maxtime_real = max_time,
            out_hints = 0,
            outlev = 0,
        ),
)
stats = bmark_solvers(solvers, cutest_problems)

using JLD2
@save "$(today())_$(prod(String.(keys(solvers)) .* "_"))cutest_$(string(length(problems)))_$(nmax).jld2" stats
