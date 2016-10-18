#export ARCSpectral, ARCSpectral_abs, ARCLDLt, ARCqKOp
export ALL_solvers
# Valid combinations
#



ALL_solvers = Symbol[] #Array(Symbol,1)


include("ARCSpectral.jl")
push!(ALL_solvers,fname)

include("ARCSpectral_abs.jl")
push!(ALL_solvers,fname)



include("TRSpectral.jl")
push!(ALL_solvers,fname)

include("TRSpectral_abs.jl")
push!(ALL_solvers,fname)



include("ARCLDLt.jl")
push!(ALL_solvers,fname)

include("ARCLDLt_abs.jl")
push!(ALL_solvers,fname)



include("TRLDLt.jl")
push!(ALL_solvers,fname)

include("TRLDLt_abs.jl")
push!(ALL_solvers,fname)



include("ARCqKOp.jl")
push!(ALL_solvers,fname)


include("ARCqKsparse.jl")
push!(ALL_solvers,fname)

include("ARCqKdense.jl")
push!(ALL_solvers,fname)



include("TRKOp.jl")
push!(ALL_solvers,fname)

include("TRKsparse.jl")
push!(ALL_solvers,fname)

include("TRKdense.jl")
push!(ALL_solvers,fname)



include("ST_TROp.jl")
push!(ALL_solvers,fname)

include("ST_TRsparse.jl")
push!(ALL_solvers,fname)

include("ST_TRdense.jl")
push!(ALL_solvers,fname)



include("ST_ARCOp.jl")
push!(ALL_solvers,fname)

include("ST_ARCsparse.jl")
push!(ALL_solvers,fname)

include("ST_ARCdense.jl")
push!(ALL_solvers,fname)



include("ARCMA97.jl")
push!(ALL_solvers,fname)

include("ARCMA97_abs.jl")
push!(ALL_solvers,fname)

include("TRMA97.jl")
push!(ALL_solvers,fname)

include("TRMA97_abs.jl")
push!(ALL_solvers,fname)

#deleteat!(ALL_solvers,1)
