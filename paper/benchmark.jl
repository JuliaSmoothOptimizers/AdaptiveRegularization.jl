using Pkg;
Pkg.activate("");
Pkg.add(url="https://github.com/JuliaSmoothOptimizers/Krylov.jl", rev="add-warm-start")
Pkg.add(url="https://github.com/vepiteski/ARCTR.jl", rev="add-callback-cg")
Pkg.update()
using CUTEst, Dates, NLPModels, SolverBenchmark

using ARCTR, JSOSolvers, NLPModelsIpopt

nmax = 1000000
problems = readlines("list_problems_$nmax.dat")
cutest_problems = (CUTEstModel(p) for p in problems)

max_time = 3600.0
max_ev = typemax(Int)
max_iter = typemax(Int)
atol = 1e-7
rtol = 1e-8

solvers = Dict(
    :ARCqKOp =>
        nlp -> ARCqKOp(
            nlp,
            verbose = false,
            atol = atol,
            rtol = rtol,
            max_time = max_time,
            max_iter = max_iter,
        ),
    :trunk =>
        nlp -> trunk(
            nlp,
            verbose = false,
            atol = atol,
            rtol = rtol,
            max_time = max_time,
            max_iter = max_iter,
        ),
)
stats = bmark_solvers(solvers, cutest_problems, skipif = prob -> (prob.meta.name in unsolved_problems))

using JLD2
@save "$(today())_$(prod(String.(keys(solvers)) .* "_"))cutest_$(string(length(problems)))_$(nmax).jld2" stats
