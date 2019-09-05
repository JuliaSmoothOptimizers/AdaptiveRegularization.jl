module ARCTR
using NLPModels

using LinearOperators
using LinearAlgebra, GenericLinearAlgebra
using Krylov
using Printf
using SparseArrays
using State
using Stopping
using Quadmath
using LDLFactorizations
using HSL

include("Includes.jl")

# include("TRARC.jl")
include("TRARCStop.jl")
# include("TRARCStop_MP.jl")
# include("TRARCStop-HO.jl")

# include("Solvers/solvers.jl")
include("Solvers/SolversStopping/solvers_stopping.jl")

end # module
