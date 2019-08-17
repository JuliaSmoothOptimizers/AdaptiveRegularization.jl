#export ARCSpectral, ARCSpectral_abs, ARCLDLt, ARCqKOp
export ALL_solvers_stopping

# Valid combinations
#

ALL_solvers_stopping = Function[]


include("ARCSpectral.jl")
push!(ALL_solvers_stopping, ARCSpectral)

include("ARCSpectral_abs.jl")
push!(ALL_solvers_stopping, ARCSpectral_abs)



include("TRSpectral.jl")
push!(ALL_solvers_stopping, TRSpectral)

include("TRSpectral_abs.jl")
push!(ALL_solvers_stopping, TRSpectral_abs)



include("ARCLDLt.jl")
push!(ALL_solvers_stopping, ARCLDLt)

include("ARCLDLt_abs.jl")
push!(ALL_solvers_stopping, ARCLDLt_abs)



include("TRLDLt.jl")
push!(ALL_solvers_stopping, TRLDLt)

include("TRLDLt_HO.jl")
push!(ALL_solvers_stopping, TRLDLt_HO)

include("TRLDLt_HO_Sham.jl")
push!(ALL_solvers_stopping, TRLDLt_HO_Sham)

include("TRLDLt_HO_Sham_BFGS.jl")
push!(ALL_solvers_stopping, TRLDLt_HO_Sham_BFGS)

include("TRLDLt_HO_MP.jl")
push!(ALL_solvers_stopping, TRLDLt_HO_MP)

include("TRLDLt_abs.jl")
push!(ALL_solvers_stopping, TRLDLt_abs)



include("ARCqKOp.jl")
push!(ALL_solvers_stopping, ARCqKOp)

include("ARCqKSparse.jl")
push!(ALL_solvers_stopping, ARCqKsparse)

include("ARCqKdense.jl")
push!(ALL_solvers_stopping, ARCqKdense)



include("TRKOp.jl")
push!(ALL_solvers_stopping, TRKOp)

include("TRKsparse.jl")
push!(ALL_solvers_stopping, TRKsparse)

include("TRKdense.jl")
push!(ALL_solvers_stopping, TRKdense)




include("ST_TROp.jl")
push!(ALL_solvers_stopping, ST_TROp)

include("ST_TROp_Sham.jl")
push!(ALL_solvers_stopping, ST_TROp_Sham)

include("ST_TROp_Sham_BFGS.jl")
push!(ALL_solvers_stopping, ST_TROp_Sham_BFGS)

include("ST_TRsparse.jl")
push!(ALL_solvers_stopping, ST_TRsparse)

include("ST_TRdense.jl")
push!(ALL_solvers_stopping, ST_TRdense)



include("ST_ARCOp.jl")
push!(ALL_solvers_stopping, ST_ARCOp)

include("ST_ARCsparse.jl")
push!(ALL_solvers_stopping, ST_ARCsparse)

include("ST_ARCdense.jl")
push!(ALL_solvers_stopping, ST_ARCdense)
#

## Will update them in the future
# include("ARCMA97.jl")
# push!(ALL_solvers,eval(fname))
#
# include("ARCMA97_abs.jl")
# push!(ALL_solvers,eval(fname))
#
# include("TRMA97.jl")
# push!(ALL_solvers,eval(fname))
#
# include("TRMA97_abs.jl")
# push!(ALL_solvers,eval(fname))
#
#
#
#
# include("ARCMA57.jl")
# push!(ALL_solvers,eval(fname))
#
# include("ARCMA57_abs.jl")
# push!(ALL_solvers,eval(fname))
#
# include("TRMA57.jl")
# push!(ALL_solvers,eval(fname))
#
# include("TRMA57_abs.jl")
# push!(ALL_solvers,eval(fname))
