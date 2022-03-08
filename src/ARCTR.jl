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

include("Solvers/solvers.jl")

end # module
