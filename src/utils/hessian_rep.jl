struct HessDense
    function HessDense(::AbstractNLPModel{T, S}, n) where {T,S}
        return new()
    end
end

struct HessSparse{S, Vi}
    rows::Vi
    cols::Vi
    vals::S
    function HessSparse(nlp::AbstractNLPModel{T, S}, n) where {T,S}
        rows, cols = hess_structure(nlp)
        vals = S(undef, nlp.meta.nnzh)
        return new{S, typeof(rows)}(rows, cols, vals)
    end
end

struct HessSparseCOO{Tv, Ti}
    H::Symmetric{Tv, SparseMatrixCOO{Tv, Ti}}
end

function HessSparseCOO(nlp::AbstractNLPModel{T, S}, n) where {T,S}
    rows, cols = hess_structure(nlp)
    vals = S(undef, nlp.meta.nnzh)
    H = Symmetric(SparseMatrixCOO(n, n, rows, cols, vals), :L)
    return HessSparseCOO(H)
end

struct HessOp{S}
    Hv::S
    function HessOp(::AbstractNLPModel{T, S}, n) where {T,S}
        return new{S}(S(undef, n))
    end
end

export HessDense, HessSparse, HessSparseCOO, HessOp

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
