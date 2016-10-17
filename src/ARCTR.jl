module ARCTR
using NLPModels
#using AmplNLReader
using OptimizationProblems
using LinearOperators
using Krylov
#using BenchmarkProfiles
using Compat
import Compat.String
# package code goes here

#using Debug

include("Includes.jl")

include("TRARC.jl")

export ARCSpectral

include("Solvers/solvers.jl")

end # module
