#export ARCSpectral, ARCSpectral_abs, ARCLDLt, ARCqKOp
export solvers
# Valid combinations
#



solvers = Array(Symbol,1)


include("ARCSpectral.jl")
push!(solvers,fname)

include("ARCSpectral_abs.jl")
push!(solvers,fname)



include("TRSpectral.jl")
push!(solvers,fname)

include("TRSpectral_abs.jl")
push!(solvers,fname)



include("ARCLDLt.jl")
push!(solvers,fname)

include("ARCLDLt_abs.jl")
push!(solvers,fname)



include("TRLDLt.jl")
push!(solvers,fname)

include("TRLDLt_abs.jl")
push!(solvers,fname)



include("ARCqKOp.jl")
push!(solvers,fname)

include("ARCqKsparse.jl")
push!(solvers,fname)

include("ARCqKdense.jl")
push!(solvers,fname)



include("TRKOp.jl")
push!(solvers,fname)

include("TRKsparse.jl")
push!(solvers,fname)

include("TRKdense.jl")
push!(solvers,fname)



include("ST_TROp.jl")
push!(solvers,fname)

include("ST_TRsparse.jl")
push!(solvers,fname)

include("ST_TRdense.jl")
push!(solvers,fname)



include("ST_ARCOp.jl")
push!(solvers,fname)

include("ST_ARCsparse.jl")
push!(solvers,fname)

include("ST_ARCdense.jl")
push!(solvers,fname)



include("ARCMA97.jl")
push!(solvers,fname)

include("ARCMA97_abs.jl")
push!(solvers,fname)

include("TRMA97.jl")
push!(solvers,fname)

include("TRMA97_abs.jl")
push!(solvers,fname)

deleteat!(solvers,1)
