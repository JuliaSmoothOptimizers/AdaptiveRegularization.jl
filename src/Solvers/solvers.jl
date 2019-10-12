export ALL_solvers

# Valid combinations

ALL_solvers = Function[]

include("ARCSpectral.jl")
push!(ALL_solvers, ARCSpectral)

include("ARCSpectral_abs.jl")
push!(ALL_solvers, ARCSpectral_abs)


include("TRSpectral.jl")
push!(ALL_solvers, TRSpectral)

include("TRSpectral_abs.jl")
push!(ALL_solvers, TRSpectral_abs)


include("ARCLDLt.jl")
push!(ALL_solvers, ARCLDLt)

include("ARCLDLt_abs.jl")
push!(ALL_solvers, ARCLDLt_abs)


include("TRLDLt.jl")
push!(ALL_solvers, TRLDLt)

include("TRLDLt_HO.jl")
push!(ALL_solvers, TRLDLt_HO)

include("TRLDLt_BBK.jl")
push!(ALL_solvers, TRLDLt_BK)

include("TRLDLt_HO_Sham.jl")
push!(ALL_solvers, TRLDLt_HO_Sham)

include("TRLDLt_HO_Sham_lambda.jl")
push!(ALL_solvers, TRLDLt_HO_Sham_λ)

include("TRLDLt_HO_Sham_BFGS.jl")
push!(ALL_solvers, TRLDLt_HO_Sham_BFGS)

# include("TRLDLt_HO_MP.jl")                # Still a work in progress
# push!(ALL_solvers, TRLDLt_HO_MP)

include("TRLDLt_abs.jl")
push!(ALL_solvers, TRLDLt_abs)


include("ARCqKOp.jl")
push!(ALL_solvers, ARCqKOp)

include("ARCqKSparse.jl")
push!(ALL_solvers, ARCqKsparse)

include("ARCqKdense.jl")
push!(ALL_solvers, ARCqKdense)


include("TRKOp.jl")
push!(ALL_solvers, TRKOp)

include("TRKsparse.jl")
push!(ALL_solvers, TRKsparse)

include("TRKdense.jl")
push!(ALL_solvers, TRKdense)


include("ST_TROp.jl")
push!(ALL_solvers, ST_TROp)

include("ST_TROp_Sham.jl")
push!(ALL_solvers, ST_TROp_Sham)

include("ST_TROp_Sham_BFGS.jl")
push!(ALL_solvers, ST_TROp_Sham_BFGS)

include("ST_TRsparse.jl")
push!(ALL_solvers, ST_TRsparse)

include("ST_TRdense.jl")
push!(ALL_solvers, ST_TRdense)


include("ST_ARCOp.jl")
push!(ALL_solvers, ST_ARCOp)

include("ST_ARCsparse.jl")
push!(ALL_solvers, ST_ARCsparse)

include("ST_ARCdense.jl")
push!(ALL_solvers, ST_ARCdense)


# include("ARCMA97.jl")
# push!(ALL_solvers, ARCMA97)
#
# include("ARCMA97_abs.jl")
# push!(ALL_solvers, ARCMA97_abs)
#
# include("TRMA97.jl")
# push!(ALL_solvers, TRMA97)
#
# include("TRMA97_abs.jl")
# push!(ALL_solvers, TRMA97_abs)

include("ARCMA57.jl")
push!(ALL_solvers, ARCMA57)

include("ARCMA57_abs.jl")
push!(ALL_solvers, ARCMA57_abs)

include("TRMA57.jl")
push!(ALL_solvers, TRMA57)

include("TRMA57_Sham.jl")
push!(ALL_solvers, TRMA57_Sham)

include("TRMA57-2.jl")
push!(ALL_solvers, TRMA57_2)

include("TRMA57_Sham_2.jl")
push!(ALL_solvers, TRMA57_Sham_2)

include("TRMA57_Sham_BFGS.jl")
push!(ALL_solvers, TRMA57_Sham_BFGS)

include("TRMA57_abs.jl")
push!(ALL_solvers, TRMA57_abs)
