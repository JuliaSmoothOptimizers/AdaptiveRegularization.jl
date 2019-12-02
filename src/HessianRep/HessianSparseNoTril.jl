export hessian_sparse_no_tril
function hessian_sparse_no_tril(nlp, x)
    # n = length(x)
    H = hess(nlp, x)
    # tempH = (H + tril(H,-1)')
    # H = tempH
    return H
end
