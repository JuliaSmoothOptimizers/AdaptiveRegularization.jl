using ARCTR
using Base.Test

# write your own tests here
@test 1 == 1
#using Optimize
using NLPModels
using AmplNLReader
using OptimizationProblems
using Compat
import Compat.String


function run_solver(solver :: Symbol, nlp :: AbstractNLPModel; args...)
  solver_f = eval(solver)
  args = Dict(args)
  skip = haskey(args, :skipif) ? pop!(args, :skipif) : x -> false
  skip(nlp) && throw(SkipException())

  # Julia nonsense
  optimal = false
  f = 0.0
  gNorm = 0.0
  status = "fail"
  #try
    (x, f, gNorm, iter, optimal, tired, status) = solver_f(nlp; args...)
  #catch e
  #  status = e.msg
  #end
  # if nlp.scale_obj
  #   f /= nlp.scale_obj_factor
  #   gNorm /= nlp.scale_obj_factor
  # end
  @printf("%-15s  %8d  %9.2e  %7.1e  %5d  %5d  %6d  %s\n",
          nlp.meta.name, nlp.meta.nvar, f, gNorm,
          nlp.counters.neval_obj, nlp.counters.neval_grad,
          nlp.counters.neval_hprod, status)
  return optimal ? (nlp.counters.neval_obj, nlp.counters.neval_grad, nlp.counters.neval_hprod) : (-nlp.counters.neval_obj, -nlp.counters.neval_grad, -nlp.counters.neval_hprod)
end



models = [AmplModel("dixmaanj.nl"), MathProgNLPModel(dixmaanj(), name="dixmaanj")]
@static if is_unix()
  using CUTEst
  push!(models, CUTEstModel("DIXMAANJ", "-param", "M=30"))
end
#solvers = [:ARCSpectral, :ARCSpectral_abs, :ARCLDLt, :ARCqKOp, :ARCqKsparse, :ARCqKdense]

for solver in solvers
  println(solver)
  for model in models
    stats = run_solver(solver, model, verbose=false)
    assert(all([stats...] .>= 0))
    reset!(model)
  end
end

# clean up the test directory
@static if is_unix()
  here = dirname(@__FILE__)
  so_files = filter(x -> (ismatch(r".so$", x) || ismatch(r".dylib$", x)), readdir(here))

  for so_file in so_files
    rm(joinpath(here, so_file))
  end

  rm(joinpath(here, "AUTOMAT.d"))
  rm(joinpath(here, "OUTSDIF.d"))
end
