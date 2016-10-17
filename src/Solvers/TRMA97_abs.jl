fname = :TRMA97_abs
c = Combi(hessian_dense,PDataMA97,solve_modelTRDiagAbs,preprocessMA97,decreaseFact,Tparam())
include("Template.jl")
