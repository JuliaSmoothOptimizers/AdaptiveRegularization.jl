using Pkg;
Pkg.activate("");
using AdaptiveRegularization, CUTEst, Dates, NLPModels, SolverBenchmark

nmax = 1000000
problems = readlines("paper/list_problems_$nmax.dat")
cutest_problems = (CUTEstModel(p) for p in problems)

unsolved_problems = [
  "CYCLOOCFLS",
  "SBRYBND",
  "BA-L49LS",
  "NONMSQRT",
  "BA-L73LS",
  "INDEF",
  "BA-L21LS",
  "FLETCBV3",
  "FLETCHBV",
  "BA-L16LS",
  "BA-L52LS",
]

max_time = 3600.0
max_ev = typemax(Int)
max_iter = typemax(Int)
atol = 1e-5
rtol = 1e-6

eL = 10.0
eU = 20.0
ψ = 10 ^ (1 // 2)
shifts = 10.0 .^ (collect(-eL:0.5:eU))

ζ = 0.5

γ₁ = 0.1
γ₂ = 5
η₁ = 0.1
η₂ = 0.75
TR = TrustRegion(
    10;
    acceptance_threshold = η₁,
    increase_threshold = η₂,
    increase_factor = γ₂,
    decrease_factor = γ₁,
)

solvers = Dict(
    :ARCqKOp =>
        nlp -> ARCqKOp(
            nlp,
            TR = TR,
            shifts = shifts,
            verbose = false,
            atol = atol,
            rtol = rtol,
            max_time = max_time,
            max_iter = max_iter,
            ζ = ζ,
        ),
    :ST_TROp =>
        nlp -> ST_TROp(
            nlp,
            TR = TR,
            verbose = false,
            atol = atol,
            rtol = rtol,
            max_time = max_time,
            max_iter = max_iter,
            ζ = ζ,
        ),
)

# Skipping unsolved problems to solve time :-)
stats = bmark_solvers(solvers, cutest_problems, skipif = prob -> (prob.meta.name in unsolved_problems))

using JLD2
@save "paper/$(today())_$(prod(String.(keys(solvers)) .* "_"))cutest_$(string(length(problems)))_$(nmax).jld2" stats
