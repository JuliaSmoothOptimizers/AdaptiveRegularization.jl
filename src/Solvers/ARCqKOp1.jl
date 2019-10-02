fname = :ARCqKOp1
shifts = 10.0.^(collect(-15.0:1.0:15.0))
c = Combi(hessian_operator,PDataK,solve_modelKARC,preprocessKARC,decreaseKARC,TparamsKARC(shifts))
include("Template.jl")
