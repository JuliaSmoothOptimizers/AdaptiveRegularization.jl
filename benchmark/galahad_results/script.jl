using Pkg
Pkg.activate("")
using DataFrames, Dates, JLD2

list_problems = "list_problems_1000000.dat"
# list_problems = "list_problems_test.dat"
problems = readlines(list_problems)

# problems = ["ARGLINA", "ARGLINB", "ARGLINC"]
npb = length(problems)

for solver in ["ARC", "TRU"]
    name_table = String[] # 
    time_table = Float64[] # Vector{Float64}(undef, npb)
    neval_obj_table = Int[] # Vector{Float64}(undef, npb)
    neval_grad_table = Int[] # Vector{Float64}(undef, npb)
    neval_hprod_table = Int[] # Vector{Float64}(undef, npb)
    objective_table = Float64[] # Vector{Float64}(undef, npb)
    dual_feas_table = Float64[] # Vector{Float64}(undef, npb)
    iter_table = Int[]

    for pb in problems
        push!(name_table, pb)
        result = readlines("$(pb)/$(solver)/" * pb * ".dat")
        neval_obj = NaN
        neval_grad = NaN
        neval_hprod = NaN
        objective = NaN
        dual_feas = NaN
        iter = NaN
        for line in result
            if occursin("Total time", line)
                push!(
                    time_table,
                    parse(Float64, replace(split(line, "=")[2], "seconds" => "")),
                )
            end
            if occursin("function evaluations", line)
                neval_obj = parse(Int, split(line, "=")[2])
            end
            if occursin("gradient evaluations", line)
                neval_grad = parse(Int, split(line, "=")[2])
            end
            if occursin("Hessian evaluations", line)
                neval_hprod = parse(Int, split(line, "=")[2])
            end
            if occursin("objective value", line)
                objective = parse(Float64, split(line, "=")[2])
            end
            if occursin("gradient norm", line)
                dual_feas = parse(Float64, split(line, "=")[2])
            end
            if occursin("major  iterations", line)
                iter = parse(Int, split(line, "=")[2])
            end
        end
        push!(neval_obj_table, neval_obj)
        push!(neval_grad_table, neval_grad)
        push!(neval_hprod_table, neval_hprod)
        push!(objective_table, objective)
        push!(dual_feas_table, dual_feas)
        push!(iter_table, iter)
    end
    df = DataFrame(
        [
            name_table,
            time_table,
            objective_table,
            iter_table,
            dual_feas_table,
            neval_obj_table,
            neval_grad_table,
            neval_hprod_table,
        ],
        [
            :name,
            :elapsed_time,
            :objective,
            :iter,
            :dual_feas,
            :neval_obj,
            :neval_grad,
            :neval_hprod,
        ],
    )

    open("$(today())_$(solver)_cutest_$(npb).dat", "w") do io
        println(io, df)
    end

    using JLD2
    @save "$(today())_$(solver)_cutest_$(npb).jld2" df
end
