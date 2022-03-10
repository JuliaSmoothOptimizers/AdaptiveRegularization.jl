# Selective includes.
include("./HessianRep/HessianDense.jl")
include("./HessianRep/HessianSparse.jl")
include("./HessianRep/HessianOp.jl")

include("./Types/Types.jl")
include("./Utilities/Utilities.jl")
include("./Utilities/ldlt_symm.jl")

path = joinpath(dirname(@__FILE__), "SolveModel")
files = filter(x -> x[(end - 2):end] == ".jl", readdir(path))
for file in files
    if file in ["krylov_aux.jl"]
        continue
    end
    include("SolveModel/" * file)
end

################################################################################
# Functions needed for high order corrections
################################################################################
pathHO = joinpath(dirname(@__FILE__), "SolveModel", "high-order-correction")
files = filter(x -> x[(end - 2):end] == ".jl", readdir(pathHO))
for file in files
    if file in []
        continue
    end
    include("SolveModel/high-order-correction/" * file)
end

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
## Model, temporary, shoudln't be here
################################################################################
# include("autodiff_high_order_model.jl")
