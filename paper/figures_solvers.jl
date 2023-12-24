using Pkg;
Pkg.activate("paper");
using JLD2, Plots, SolverBenchmark, DataFrames

name = "2022-05-16_ST_TROp_ARCqKOpShift05_cutest_277_1000000"
@load "paper/$name.jld2" stats

solved(df) = (df.status .== :first_order) .| (df.status .== :unbounded)

for solver in keys(stats)
  open("paper/$(name)_result_$(solver).dat", "w") do io
    print(
      io,
      stats[solver][
        !,
        [
          :name,
          :nvar,
          # :ncon,
          :status,
          :objective,
          :elapsed_time,
          :iter,
          # :primal_feas,
          :dual_feas,
          :neval_obj,
          :neval_grad,
          :neval_hprod,
          :neval_hess,
        ],
      ],
    )
  end
end

nmins = [0, 100, 1000, 10000]
for nmin in nmins
  # Same figure with minimum number of variables
  stats2 = copy(stats)
  for solver in keys(stats)
    stats2[solver] = stats[solver][stats[solver].nvar .>= nmin, :]
  end

  nb_problems = length(stats2[first(keys(stats))][!, :name])

  # Figures comparing two results:
  costs_all = [
    df -> .!solved(df) * Inf + df.elapsed_time,
    df -> .!solved(df) * Inf + df.neval_obj,
    df -> .!solved(df) * Inf + df.neval_grad,
    df -> .!solved(df) * Inf + df.neval_hprod,
  ]
  costnames_all = ["elapsed time", "objective evals", "gradient evals", "hessian-vector products"]
  p =
    profile_solvers(stats2, costs_all, costnames_all, height = 400, width = 400, margin = 5Plots.mm)
  png(p, "paper/$(name)_all($(nb_problems))_min_$(nmin)")
end
