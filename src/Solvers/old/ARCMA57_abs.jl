fname = :ARCMA57_abs
c = Combi(hessian_sparse,PDataMA57,solve_modelARCDiagAbs,preprocessMA57,decreaseFact,Tparam())
include("Template.jl")
