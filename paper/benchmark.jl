using Pkg;
Pkg.activate("");
Pkg.add(url="https://github.com/vepiteski/ARCTR.jl", rev="main")
Pkg.update()
using CUTEst, Dates, NLPModels, SolverBenchmark

using ARCTR, JSOSolvers, NLPModelsIpopt

nmax = 1000000
problems = readlines("list_problems_$nmax.dat")
cutest_problems = (CUTEstModel(p) for p in problems)

max_time = 1200.0
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
            shifts = 10.0 .^ (collect(-10.0:0.5:20.0)),
        ),
    :trunk =>
        nlp -> trunk(
            nlp,
            verbose = 0,
            atol = atol,
            rtol = rtol,
            max_time = max_time,
        ),
    :ipopt => 
        nlp -> ipopt(
            nlp,
            print_level = 0,
            dual_inf_tol = Inf,
            constr_viol_tol = Inf,
            compl_inf_tol = Inf,
            acceptable_iter = 0,
            max_cpu_time = max_time,
            tol = rtol,
        ),
)
stats = bmark_solvers(solvers, cutest_problems)

using JLD2
@save "$(today())_$(prod(String.(keys(solvers)) .* "_"))cutest_$(string(length(problems)))_$(nmax).jld2" stats
