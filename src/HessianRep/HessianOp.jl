function hessian_operator(nlp,x)
    return LinearOperator(length(x), Float64, v -> hprod(nlp,x,v))
end
