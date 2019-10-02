fname = :TRLDLt
c = Combi(hessian_dense,PDataLDLt,solve_modelTRDiag,preprocessLDLt,decreaseFact,Tparam())
include("Template.jl")
