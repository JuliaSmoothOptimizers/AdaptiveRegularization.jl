using Pkg;
Pkg.activate(".");
using JLD2, Plots, SolverBenchmark

name = "2022-06-24_trunk_ARCqKOp_ipopt_cutest_277_1000000"
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

nmin = 100
stats2 = copy(stats)
for solver in keys(stats)
    stats2[solver] = stats[solver][stats[solver].nvar.>=nmin, :]
end
pop!(stats2, :trunk)

# Number of problems solved by ipopt
@show sum(solved(stats2[:ipopt])) # 97
# Number of problems solved by dci
@show sum(solved(stats2[:ARCqKOp])) # 101

tim_dci = stats2[:ARCqKOp][solved(stats2[:ARCqKOp]), :elapsed_time]
tim_ipopt = stats2[:ipopt][solved(stats2[:ARCqKOp]), :elapsed_time]
# Number of problems where Ipopt is fastest
@show sum(tim_dci .> tim_ipopt) # 43
# Number of problems where DCI is fastest
@show sum(tim_dci .< tim_ipopt) # 58

obj_dci = stats2[:ARCqKOp][solved(stats2[:ARCqKOp]), :neval_obj]
obj_ipopt = stats2[:ipopt][solved(stats2[:ARCqKOp]), :neval_obj]
con_dci = stats2[:ARCqKOp][solved(stats2[:ARCqKOp]), :neval_grad]
con_ipopt = stats2[:ipopt][solved(stats2[:ARCqKOp]), :neval_grad]
# Number of problems where Ipopt use less evaluations
@show sum((con_dci) .> (con_ipopt)) # 40
# Number of problems where DCI use less evaluations
@show sum((con_dci) .< (con_ipopt)) # 51
# Number where it is a tie
@show sum((con_dci) .== (con_ipopt)) # 10

costs =
    [df -> .!solved(df) * Inf + df.elapsed_time, df -> .!solved(df) * Inf + df.neval_grad]
costnames = ["Time", "Number of gradient evaluations"]
p = profile_solvers(stats2, costs, costnames)
png(p, "$name")

for solver in keys(stats2)
    open("$(name)_result_$(solver).dat", "w") do io
        print(
            io,
            stats2[solver][
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
