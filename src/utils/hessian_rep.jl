"""
    HessDense(::AbstractNLPModel{T,S}, n)
Return a structure used for the evaluation of dense Hessian matrix.
"""
struct HessDense
    function HessDense(::AbstractNLPModel{T,S}, n) where {T,S}
        return new()
    end
end

"""
    HessSparse(::AbstractNLPModel{T,S}, n)
Return a structure used for the evaluation of sparse Hessian matrix.
"""
struct HessSparse{S,Vi}
    rows::Vi
    cols::Vi
    vals::S
    function HessSparse(nlp::AbstractNLPModel{T,S}, n) where {T,S}
        rows, cols = hess_structure(nlp)
        vals = S(undef, nlp.meta.nnzh)
        return new{S,typeof(rows)}(rows, cols, vals)
    end
end

"""
    HessSparseCOO(::AbstractNLPModel{T,S}, n)
Return a structure used for the evaluation of sparse Hessian matrix in COO-format.
"""
struct HessSparseCOO{Tv,Ti}
    H::Symmetric{Tv,SparseMatrixCOO{Tv,Ti}}
end

function HessSparseCOO(nlp::AbstractNLPModel{T,S}, n) where {T,S}
    rows, cols = hess_structure(nlp)
    vals = S(undef, nlp.meta.nnzh)
    H = Symmetric(SparseMatrixCOO(n, n, rows, cols, vals), :L)
    return HessSparseCOO(H)
end

"""
    HessOp(::AbstractNLPModel{T,S}, n)
Return a structure used for the evaluation of the Hessian matrix as an operator.
"""
struct HessOp{S}
    Hv::S
    function HessOp(::AbstractNLPModel{T,S}, n) where {T,S}
        return new{S}(S(undef, n))
    end
end

export HessDense, HessSparse, HessSparseCOO, HessOp

"""
    hessian!(workspace::HessDense, nlp, x)

Return the Hessian matrix of `nlp` at `x` in-place with memory update of `workspace`.
"""
function hessian! end

function hessian!(workspace::HessDense, nlp, x)
    H = Matrix(hess(nlp, x))
    return H
end

function hessian!(workspace::HessOp, nlp, x)
    return hess_op!(nlp, x, workspace.Hv)
end

function hessian!(workspace::HessSparse, nlp, x)
    hess_coord!(nlp, x, workspace.vals)
    n = nlp.meta.nvar
    return Symmetric(sparse(workspace.rows, workspace.cols, workspace.vals, n, n), :L)
end

function hessian!(workspace::HessSparseCOO, nlp, x)
    hess_coord!(nlp, x, workspace.H.data.vals)
    return workspace.H
end
