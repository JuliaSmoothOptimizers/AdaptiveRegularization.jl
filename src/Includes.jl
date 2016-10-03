# Selective includes.
include("./HessianRep/HessianDense.jl")
include("./HessianRep/HessianSparse.jl")
include("./HessianRep/HessianOp.jl")
using HSL
using LinearOperators

include("./Types/Types.jl")
include("./Utilities/Utilities.jl")

include("./SolveModel/SolveModelKARC.jl")
include("./SolveModel/SolveModelKTR.jl")
include("./SolveModel/SolveModelARCDiagAbs.jl")
include("./SolveModel/SolveModelARCDiag.jl")
include("./SolveModel/SolveModelTRDiag.jl")
include("./SolveModel/SolveModelTRDiagAbs.jl")
include("./SolveModel/SolveDiag.jl")
include("./SolveModel/SolveDiagTR.jl")
include("./SolveModel/SolveModelST_TR.jl")
include("./SolveModel/SolveModelST_ARC.jl")
include("./SolveModel/cgARC.jl")
include("./SolveModel/krylov_aux.jl")

include("./Utilities/ldlt_symm.jl")
include("./PreProcess/TParamsKARC.jl")
include("./PreProcess/TParamsKTR.jl")
include("./PreProcess/TParamsST.jl")
include("./PreProcess/PreProcessLDLt.jl")
include("./PreProcess/PreProcessSpectral.jl")
include("./PreProcess/PreProcessMA97.jl")
include("./PreProcess/PreProcessKARC.jl")
include("./PreProcess/PreProcessKTR.jl")
include("./PreProcess/PreProcessST_TR.jl")

