module ARCTR
using NLPModels

using LinearOperators
using LinearAlgebra, GenericLinearAlgebra
using Krylov
using Printf
using SparseArrays
using State
using Stopping


include("Includes.jl")

# include("TRARC.jl")
include("TRARCStop.jl")
include("TRARCStop-HO.jl")
include("TRARCStop-HO-4.jl")

# include("Solvers/solvers.jl")
include("Solvers/SolversStopping/solvers_stopping.jl")

end # module
