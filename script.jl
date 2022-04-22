using Pkg
Pkg.activate(".")
using CUTEst, NLPModels
using ARCTR # https://github.com/tmigot/ARCTR.jl#expand-param
# The problematic function:
# https://github.com/tmigot/ARCTR.jl/blob/ef38375e49735a5440090562a10d553516c231e4/src/PreProcess/PreProcessKARC.jl#L38

max_time = 600.0 #20 minutes
max_ev = typemax(Int)
max_iter = 15000
tol = 1e-5

ARCqKOpT = nlp -> ARCqKOp(
    nlp,
    shifts = 10.0 .^ (collect(-20.0:1.0:20.0)),
    verbose = true,
    atol = tol,
    rtol = tol,
    max_time = max_time,
    max_iter = max_iter,
)

nlp = CUTEstModel("NELSONLS")
stats = ARCqKOpT(nlp)
print(stats)

finalize(nlp)