fname = :TRMA57
c = Combi(hessian_sparse,PDataMA57,solve_modelTRDiag,preprocessMA57,decreaseFact,Tparam())
include("Template.jl")
