fname = :ST_TRsparse
c = Combi(hessian_sparse,PDataST,solve_modelST_TR,preprocessST,decreaseGen,TparamsST())
include("Template.jl")
