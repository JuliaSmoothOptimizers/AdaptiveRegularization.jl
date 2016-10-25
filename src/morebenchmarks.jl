# Benchmark two solvers on a set of small problems and profile.
using Optimize

# Common setup
n = 5000
mpb_probs = filter(name -> name != :OptimizationProblems && name != :sbrybnd, names(OptimizationProblems))
ampl_prob_dir = "/home/local/USHERBROOKE/dusj1701/Documents/Recherche/TR-ARC/UnifiedImplementation/Julia/AMPLUnconstrained/"
ampl_probs = [Symbol(split(p, ".")[1]) for p in filter(x -> contains(x, ".nl"), readdir(ampl_prob_dir))]

# Example 1: benchmark two solvers on a set of problems
function two_solvers()
  solvers = [trunk, lbfgs]
  title = @sprintf("f+g+hprod on %d problems of size about %d", length(mpb_probs), n)
  bmark_args = Dict{Symbol, Any}(:format => :mpb, :skipif => model -> model.meta.ncon > 0)
  profile_args = Dict{Symbol, Any}(:title => title)
  profiles = bmark_and_profile(solvers, mpb_probs, n, bmark_args=bmark_args, profile_args=profile_args)
end

# Example 1b: benchmark solvers on a set of problems
#solvers = [ARCLDLt_abs, trunk, lbfgs]
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
    #, :Jac, :Hess, :HessVec, :ExprGraph]),

  return stats, profiles
end


# Example 2: benchmark one solver on problems written in two modeling languages
function two_languages()
  probs = ampl_probs ∩ mpb_probs
  title = @sprintf("f+g+hprod on %d problems of size about %d", length(probs), n)
  stats = Dict{Symbol, Array{Int,2}}()
  stats[:trunk_ampl] = run_problems(:trunk, probs, n, format=:ampl, skipif=model -> model.meta.ncon != 0)
  stats[:trunk_mpb] = run_problems(:trunk, probs, n, format=:mpb, skipif=model -> model.meta.ncon != 0)
  profile_solvers(stats, title=title)
end

# Example 3: benchmark one solver with different options on a set of problems
function solver_options()
  title = @sprintf("f+g+hprod on %d problems of size about %d", length(mpb_probs), n)
  stats = Dict{Symbol, Array{Int,2}}()
  stats[:trunk] = run_problems(:trunk, mpb_probs, n, format=:mpb, skipif=model -> model.meta.ncon != 0, nm_itmax=5)
  stats[:trunk_monotone] = run_problems(:trunk, mpb_probs, n, format=:mpb, skipif=model -> model.meta.ncon != 0, monotone=true)
  profile_solvers(stats, title=title)
end
