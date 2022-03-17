module ARCTR

# stdlib
using LinearAlgebra, Logging, Printf, SparseArrays
# JSO
using Krylov, LinearOperators, NLPModels, SolverCore, SolverTools
# Stopping
using Stopping, StoppingInterface

# Selective includes.
include("hessian_rep.jl")

include("./Utilities/show.jl")
include("./Utilities/pdata_struct.jl")
include("./Utilities/Utilities.jl")
include("./Utilities/increase_decrease.jl")
include("./Utilities/ldlt_symm.jl")

path = joinpath(dirname(@__FILE__), "SolveModel")
files = filter(x -> x[(end-2):end] == ".jl", readdir(path))
for file in files
    if file in ["krylov_aux.jl"]
        continue
    end
    include("SolveModel/" * file)
end

pathHO = joinpath(dirname(@__FILE__), "SolveModel", "high-order-correction")
files = filter(x -> x[(end-2):end] == ".jl", readdir(pathHO))
for file in files
    if file in []
        continue
    end
    include("SolveModel/high-order-correction/" * file)
end

path = joinpath(dirname(@__FILE__), "PreProcess")
files = filter(x -> x[(end-2):end] == ".jl", readdir(path))
for file in files
    if occursin(r"MA", file)
        continue # Remove MA57 and MA97 solvers for now
    end
    include("PreProcess/" * file)
end

####################################################################
## Model, temporary, shoudln't be here
# include("autodiff_high_order_model.jl")
####################################################################

include("TRARC.jl")

export ALL_solvers

ALL_solvers = Function[]

path = joinpath(dirname(@__FILE__), "Solvers")
files = filter(x -> x[(end-2):end] == ".jl", readdir(path))
for file in files
    if occursin(r"MA", file)
        continue # Remove MA57 and MA97 solvers for now
    end
    if file in ["ARCqKOp05.jl", "ARCqKOp2.jl", "TRLDLt_HO_MP.jl"] #use something else than TRARC
        continue
    end
    include("Solvers/" * file)
    fun = Symbol(split(file, ".")[1])
    push!(ALL_solvers, eval(fun))

    @eval begin
        function $fun(nlp::AbstractNLPModel; kwargs...)
            nlpstop = NLPStopping(nlp; kwargs...)
            nlpstop = $fun(nlpstop; kwargs...)
            return stopping_to_stats(nlpstop)
        end
    end
    @eval export $fun
end

end # module
