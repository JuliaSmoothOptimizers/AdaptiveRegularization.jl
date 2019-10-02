fname = :ST_ARCsparse
c = Combi(hessian_sparse,PDataST,solve_modelST_ARC,preprocessST,decreaseGen,TparamsST())
include("Template.jl")
