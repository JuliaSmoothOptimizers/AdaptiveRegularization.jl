# Benchmark two solvers on a set of small problems and profile.
using Optimize

# Common setup
n = 5000
mpb_probs = filter(name -> name != :OptimizationProblems && name != :sbrybnd, names(OptimizationProblems))
ampl_prob_dir = "/home/local/USHERBROOKE/dusj1701/Documents/Recherche/TR-ARC/UnifiedImplementation/Julia/AMPLUnconstrained/"
ampl_probs = [Symbol(split(p, ".")[1]) for p in filter(x -> contains(x, ".nl"), readdir(ampl_prob_dir))]


# Example 1b: benchmark solvers on a set of problems
function compare_solvers(solvers,probs;title:: String = " ", kwargs...)
  bmark_args = Dict{Symbol, Any}(:skipif => model -> model.meta.ncon > 0, :max_f => 20000)#,:atol=>1e-10,:rtol=>1e-8)
  profile_args = Dict{Symbol, Any}(:title => title)
  stats, profiles = bmark_and_profile(solvers,
                                      #(MathProgNLPModel(eval(p)(n), name=string(p)) for p in mpb_probs),

                                      (MathProgNLPModel(eval(p)(n), 
                                                        name=string(p)#, 
                                                        #features = [:Grad]
                                                        ) 
                                       for p in probs),
                                      bmark_args=bmark_args, profile_args=profile_args)
    # Other available features :Jac, :Hess, :HessVec, :ExprGraph

  return stats, profiles
end
