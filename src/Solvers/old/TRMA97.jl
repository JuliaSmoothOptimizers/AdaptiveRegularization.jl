fname = :TRMA97
c = Combi(hessian_dense,PDataMA97,solve_modelTRDiag,preprocessMA97,decreaseFact,Tparam())
include("Template.jl")
