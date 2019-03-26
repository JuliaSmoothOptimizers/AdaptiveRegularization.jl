function hessian_operator(nlp, x)
    n = nlp.meta.nvar
    # temp = Array(Float64, n)
    temp = Array{Float64}(undef, n)
    return hess_op!(nlp, x, temp)
    #return LinearOperator(nlp.meta.nvar, nlp.meta.nvar, true, true, v -> hprod(nlp,x,v))
end
