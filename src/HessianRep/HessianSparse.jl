function hessian_sparse(nlp,x)
    n = length(x)
    H=hess(nlp,x)
    tempH = (H+tril(H,-1)')
    H = tempH
    return H
end

