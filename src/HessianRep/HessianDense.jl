export hessian_dense
function hessian_dense(nlp,x)
    n = length(x)
    H = hess(nlp, x)
    tempH = (H + tril(H, -1)')
    H = Matrix(Symmetric(tempH, :L))
    return H
end
