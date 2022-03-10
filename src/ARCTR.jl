module ARCTR

using NLPModels
using LinearOperators
using LinearAlgebra
using Krylov
using Printf
using SparseArrays
using Stopping
using Logging, SolverCore, SolverTools

include("Includes.jl")

include("TRARC.jl")

####################################################################
# Include solvers
export ALL_solvers

ALL_solvers = Function[]

path = joinpath(dirname(@__FILE__), "Solvers")
files = filter(x -> x[(end - 2):end] == ".jl", readdir(path))
for file in files
    if occursin(r"MA", file)
        continue # Remove MA57 and MA97 solvers for now
    end
    if file in ["ARCqKOp05.jl", "ARCqKOp2.jl", "TRLDLt_HO_MP.jl"] #use something else than TRARC
        continue
    end
    include("Solvers/" * file)
    push!(ALL_solvers, eval(Symbol(split(file, ".")[1])))
end
####################################################################

end # module
