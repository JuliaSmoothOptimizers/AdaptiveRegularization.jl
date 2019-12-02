fname = :ARCqKsparse
c = Combi(hessian_sparse,PDataK,solve_modelKARC,preprocessKARC,decreaseKARC,TparamsKARC())
include("Template.jl")
