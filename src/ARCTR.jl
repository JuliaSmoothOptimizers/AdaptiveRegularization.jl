module ARCTR
using NLPModels

using LinearOperators
using LinearAlgebra
using Krylov
using Printf
using SparseArrays
using State
using Stopping


include("Includes.jl")

include("TRARC.jl")
include("TRARCStop2.jl")

include("Solvers/solvers.jl")
include("Solvers/SolversStopping/solvers_stopping.jl")

end # module
