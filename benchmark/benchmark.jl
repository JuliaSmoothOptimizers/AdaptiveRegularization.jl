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
  :ipopt => nlp -> ipopt(
    nlp,
    print_level = 0,
    dual_inf_tol = Inf,
    constr_viol_tol = Inf,
    compl_inf_tol = Inf,
    acceptable_iter = 0,
    max_cpu_time = max_time,
    tol = rtol,
  ),
  :lbfgs => nlp -> lbfgs(
    nlp,
    verbose = 0,
    atol = atol,
    rtol = rtol,
    max_time = max_time,
    max_eval = max_ev,
  ),
  :trunk => nlp -> trunk(
    nlp,
    verbose = 0,
    atol = atol,
    rtol = rtol,
    max_time = max_time,
    max_eval = max_ev,
  ),
  :tron => nlp -> tron(
    nlp,
    verbose = 0,
    atol = atol,
    rtol = rtol,
    max_time = max_time,
    max_eval = max_ev,
  ),
)
stats = bmark_solvers(solvers, cutest_problems)

using JLD2
@save "$(today())_$(prod(String.(keys(solvers)) .* "_"))cutest_$(string(length(problems)))_$(nmax).jld2" stats
