export hessian_dense
function hessian_dense(nlp,x)
    H = Matrix(hess(nlp, x))
    return H
end
