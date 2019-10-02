fname = :TRKOp
c = Combi(hessian_operator, PDataK, solve_modelKTR, preprocessKTR, decreaseKTR, TparamsKTR())
include("Template.jl")
