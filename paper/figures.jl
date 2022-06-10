using Pkg;
Pkg.activate(".");
using JLD2, Plots, SolverBenchmark

name = "..."
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
    [df -> .!solved(df) * Inf + df.elapsed_time, df -> .!solved(df) * Inf + df.neval_hprod]
costnames = ["Time", "Number of Hessian-vector products"]
p = profile_solvers(stats, costs, costnames)
png(p, "$name")

for solver in keys(stats)
    open("$(name)_result_$(solver).dat", "w") do io
        print(
            io,
            stats[solver][
                !,
                [
                    :name,
                    :nvar,
                    :ncon,
                    :status,
                    :objective,
                    :elapsed_time,
                    :iter,
                    :primal_feas,
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
