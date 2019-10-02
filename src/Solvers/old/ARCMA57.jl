fname = :ARCMA57
c = Combi(hessian_sparse,PDataMA57,solve_modelARCDiag,preprocessMA57,decreaseFact,Tparam())
include("Template.jl")
