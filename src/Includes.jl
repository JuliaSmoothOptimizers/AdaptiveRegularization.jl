# Selective includes.
include("./HessianRep/HessianDense.jl")
include("./HessianRep/HessianSparse.jl")
include("./HessianRep/HessianSparseNoTril.jl")
include("./HessianRep/HessianOp.jl")

include("./Types/Types.jl")
include("./Utilities/Utilities.jl")

include("./SolveModel/SolveModelKARC.jl")
include("./SolveModel/SolveModelKTR.jl")
include("./SolveModel/SolveModelARCDiagAbs.jl")
include("./SolveModel/SolveModelARCDiag.jl")
include("./SolveModel/SolveModelARCDiag_HO.jl")
include("./SolveModel/SolveModelARCDiag_HO_vs_Nwt.jl")
include("./SolveModel/SolveModelTRDiag.jl")
include("./SolveModel/SolveModelTRDiagAbs.jl")
include("./SolveModel/SolveModelTRDiag-HO-lambda-vs-Nwt.jl")
include("./SolveModel/SolveDiag.jl")
include("./SolveModel/SolveDiagTR.jl")
include("./SolveModel/SolveModelST_TR.jl")
include("./SolveModel/SolveModelST_ARC.jl")
include("./SolveModel/cgARC.jl")
#include("./SolveModel/krylov_aux.jl")

include("./Utilities/ldlt_symm.jl")
include("./PreProcess/TParamsKARC.jl")
include("./PreProcess/TParamsKTR.jl")
include("./PreProcess/TParamsST.jl")
include("./PreProcess/PreProcessLDLt.jl")
include("./PreProcess/PreProcessSpectral.jl")
# include("./PreProcess/PreProcessMA97.jl")
# include("./PreProcess/PreProcessMA57.jl")
include("./PreProcess/PreProcessKARC.jl")
include("./PreProcess/PreProcessKTR.jl")
include("./PreProcess/PreProcessST_TR.jl")


################################################################################
# Functions needed for high order corrections
################################################################################
# include("./SolveModel/SolveModelTRDiag-HO.jl")
# include("./SolveModel/SolveModelTRDiag-HO-lambda.jl")
include("./SolveModel/high-order-correction/shamanskii.jl")
include("./SolveModel/high-order-correction/shamanskii-lambda.jl")
# include("./SolveModel/high-order-correction/shamanskii-ma57.jl")
# include("./SolveModel/high-order-correction/shamanskii-ma57-bfgs.jl")
include("./SolveModel/high-order-correction/shamanskii-bfgs.jl")
include("./SolveModel/high-order-correction/chebyshev.jl")
include("./SolveModel/high-order-correction/halley.jl")
include("./SolveModel/high-order-correction/super-halley.jl")


################################################################################
## Model, temporary, shoudln't be here
################################################################################
# include("autodiff_high_order_model.jl")
