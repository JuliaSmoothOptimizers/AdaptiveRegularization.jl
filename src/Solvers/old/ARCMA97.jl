fname = :ARCMA97
c = Combi(hessian_dense,PDataMA97,solve_modelARCDiag,preprocessMA97,decreaseFact,Tparam())
include("Template.jl")
