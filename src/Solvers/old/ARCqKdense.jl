fname = :ARCqKdense
c = Combi(hessian_dense,PDataK,solve_modelKARC,preprocessKARC,decreaseKARC,TparamsKARC())
include("Template.jl")
