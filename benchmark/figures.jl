using Pkg;
Pkg.activate(".");
using JLD2, Plots, SolverBenchmark

name = "2022-03-20_ST_TROp_ARCqKOp_cutest_277_1000000"
@load "$name.jld2" stats
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

if collect(keys(stats)) == [:ST_TROp, :ARCqKOp]
    # Same figure with minimum number of variables
    nmin = 100
    stats2 = copy(stats)
    stats2[:ST_TROp] = stats[:ST_TROp][stats[:ST_TROp].nvar.>nmin, :]
    stats2[:ARCqKOp] = stats[:ARCqKOp][stats[:ARCqKOp].nvar.>nmin, :]

    p = profile_solvers(stats2, costs, costnames)
    png(p, "$(name)_min_$(nmin + 1)")

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
    png(p, "$(name)_all_min_$(nmin + 1)")
    p = profile_solvers(stats, costs_all, costnames_all)
    png(p, "$(name)_all")

end
