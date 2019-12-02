export hessian_operator
function hessian_operator(nlp, x)
    n = nlp.meta.nvar
    temp = Array{Float64}(undef, n)
    return hess_op!(nlp, x, temp)
end
