using ARCTR
using Base.Test
using NLPModels
using JuMP

# test with the well known Woods test function
include("woods.jl")

function run_solver(solver_f :: Function, nlp :: AbstractNLPModel; args...)
  #solver_f = eval(solver)
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


#########################################  tests begin here #########################

include("../src/Utilities/testLDLt.jl")

models = [MathProgNLPModel(woods(), name="woods")]

nbsolver = 0
for solver in ALL_solvers
  nbsolver += 1
  println(nbsolver,"  ",solver)
  for model in models
    stats = run_solver(solver, model, verbose=false)
    @test (all([stats...] .>= 0))
    reset!(model)
  end
end
