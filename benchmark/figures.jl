using Pkg;
Pkg.activate(".");
using JLD2, Plots, SolverBenchmark

name = "2022-04-04_ARCqKOp10201z075_ARCqKOp10201z05_ARCqKOp10201z025_ARCqKOp10201z1_cutest_263_10000"
@load "$name.jld2" stats
#=
stats1 = copy(stats)
name = "2022-03-21_trunk_tron_ipopt_lbfgs_cutest_277_1000000"
@load "$name.jld2" stats
stats = merge(stats, stats1)
=#
solved(df) = (df.status .== :first_order) .| (df.status .== :unbounded)

open("$name.dat", "w") do io
    for solver in keys(stats)
        # Number of problems solved
        println(io, solver)
        println(io, size(stats[solver][solved(stats[solver]), [:name]], 1))
        println(
            io,
            stats[solver][
                .!solved(stats[solver]),
                [:name, :nvar, :status, :elapsed_time, :dual_feas],
            ],
        )
    end
end

costs =
    [df -> .!solved(df) * Inf + df.elapsed_time, df -> .!solved(df) * Inf + df.neval_obj]
costnames = ["Time", "Evaluations of f"]
p = profile_solvers(stats, costs, costnames)
png(p, "$name")

#=
open("$name.dat", "w") do io
  print(io, stats[:fps][!, [:name, :nvar, :ncon, :status, :objective, :elapsed_time, :iter, :primal_feas, :dual_feas]])
end
=#

nmins = [0, 100, 1000, 10000]
for nmin in nmins
    # Same figure with minimum number of variables
    stats2 = copy(stats)
    for solver in keys(stats)
        stats2[solver] = stats[solver][stats[solver].nvar.>=nmin, :]
    end

    nb_problems = length(stats2[first(keys(stats))][!, :name])

    # Figures comparing two results:
    costs_all = [
        df -> .!solved(df) * Inf + df.elapsed_time,
        df -> .!solved(df) * Inf + df.neval_obj,
        df -> .!solved(df) * Inf + df.neval_grad,
        df -> .!solved(df) * Inf + df.neval_hprod,
        df -> .!solved(df) * Inf + df.neval_obj + df.neval_grad + df.neval_hprod,
    ]
    costnames_all = [
        "time",
        "objective evals",
        "gradient evals",
        "hessian-vector products",
        "obj + grad + hprod",
    ]
    p = profile_solvers(stats2, costs_all, costnames_all)
    png(p, "$(name)_all($(nb_problems))_min_$(nmin)")
end
