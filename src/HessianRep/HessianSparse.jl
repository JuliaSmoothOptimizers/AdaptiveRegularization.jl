export hessian_sparse
function hessian_sparse(nlp, x)
    H = hess(nlp, x)
    return H
end
