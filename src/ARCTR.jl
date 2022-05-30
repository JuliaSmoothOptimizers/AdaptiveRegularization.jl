module ARCTR

# stdlib
using LinearAlgebra, Logging, Printf, SparseArrays
# JSO
using HSL, Krylov, LinearOperators, NLPModels, SolverCore, SolverTools
# Stopping
using Stopping, StoppingInterface

# Selective includes.
include("./utils/hessian_rep.jl")
include("./utils/pdata_struct.jl")
include("./utils/utils.jl")
include("./utils/increase_decrease.jl")

path = joinpath(dirname(@__FILE__), "SolveModel")
files = filter(x -> x[(end-2):end] == ".jl", readdir(path))
for file in files
    if file in ["krylov_aux.jl"]
        continue
    end
    include("SolveModel/" * file)
end

path = joinpath(dirname(@__FILE__), "PreProcess")
files = filter(x -> x[(end-2):end] == ".jl", readdir(path))
for file in files
    if occursin(r"MA", file)
        continue # Remove MA57 and MA97 solvers for now
    end
    include("PreProcess/" * file)
end

include("main.jl")

export ALL_solvers

ALL_solvers = Symbol[]

include("solvers.jl")

for fun in keys(solvers_const)
    push!(ALL_solvers, fun)

    ht, pt, sm, ka = ARCTR.solvers_const[fun]
    @eval begin
        function $fun(nlpstop::NLPStopping; kwargs...)
            kw_list = Dict{Symbol,Any}()
            if $ka != ()
                for t in $ka
                    push!(kw_list, t)
                end
            end
            merge!(kw_list, Dict(kwargs))
            TRARC(nlpstop; hess_type = $ht, pdata_type = $pt, solve_model = $sm, kw_list...)
        end
    end
    @eval begin
        function $fun(nlp::AbstractNLPModel; kwargs...)
            # https://github.com/vepiteski/Stopping.jl/blob/4036e996be4e9521bb40725a659a03d23fd70e03/src/Stopping/nlp_admissible_functions.jl#L12
            nlpstop = NLPStopping(
                nlp;
                optimality_check = (pb, state) ->
                    unconstrained_check(pb, state, pnorm = 2.0),
                kwargs...,
            )
            nlpstop = $fun(nlpstop; kwargs...)
            return stopping_to_stats(nlpstop)
        end
    end
    @eval export $fun
end

end # module
