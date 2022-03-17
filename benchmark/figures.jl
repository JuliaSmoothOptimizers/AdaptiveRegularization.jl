using Pkg; Pkg.activate(".")
using JLD2, Plots, SolverBenchmark

name = "2022-03-17_ST_TROp_ARCqKOp_cutest_129_10"
@load "$name.jld2" stats
solved(df) = (df.status .== :first_order) .| (df.status .== :unbounded)

open("$name.dat", "w") do io
    for solver in keys(stats)
    # Number of problems solved
    println(io, solver)
    println(io, size(stats[solver][solved(stats[solver]), [:name]], 1))
    println(io, stats[solver][.!solved(stats[solver]), [:name, :nvar, :status, :elapsed_time, :dual_feas]])
    end
end

costs = [
  df -> .!solved(df) * Inf + df.elapsed_time,
  df -> .!solved(df) * Inf + df.neval_obj,
]
costnames = ["Time", "Evaluations of f"]
p = profile_solvers(stats, costs, costnames)
png(p, "$name")

#=
open("$name.dat", "w") do io
  print(io, stats[:fps][!, [:name, :nvar, :ncon, :status, :objective, :elapsed_time, :iter, :primal_feas, :dual_feas]])
end
=#
