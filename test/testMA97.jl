using ARCTR
using Base.Test
using NLPModels
using JuMP

# first test the LDLt function
include("../src/Utilities/testLDLt.jl")

# test all solvers with the well known Woods test function
include("woods.jl")
nlp = MathProgNLPModel(woods(), name="woods")

solver = ALL_solvers[21]
(x, f, gNorm, iter, optimal, tired, status) = solver(nlp, verbose=true)

reset!(nlp);
