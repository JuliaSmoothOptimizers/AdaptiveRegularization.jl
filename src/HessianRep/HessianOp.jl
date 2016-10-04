function hessian_operator(nlp,x)
    return LinearOperator(nlp.meta.nvar, nlp.meta.nvar, true, true, v -> hprod(nlp,x,v))
end
