abstract type AbstractHess end

function get_hess(Hstruct::AbstractHess)
  Hstruct.H
end

"""
    HessDense(::AbstractNLPModel{T,S}, n)
Return a structure used for the evaluation of dense Hessian matrix.
"""
struct HessDense{T} <: AbstractHess
  H::Matrix{T}
  function HessDense(::AbstractNLPModel{T, S}, n) where {T, S}
    H = Matrix{Float64}(undef, n, n)
    return new{T}(H)
  end
end

"""
    HessSparse(::AbstractNLPModel{T,S}, n)
Return a structure used for the evaluation of sparse Hessian matrix.
"""
struct HessSparse{T, S, Vi, It <: Integer} <: AbstractHess
  rows::Vi
  cols::Vi
  vals::S
  H::Symmetric{T, SparseMatrixCSC{T, It}}
  function HessSparse(nlp::AbstractNLPModel{T, S}, n) where {T, S}
    rows, cols = hess_structure(nlp)
    vals = S(undef, nlp.meta.nnzh)
    H = Symmetric(spzeros(T, n, n), :L)
    return new{T, S, typeof(rows), eltype(rows)}(rows, cols, vals, H)
  end
end

"""
    HessSparseCOO(::AbstractNLPModel{T,S}, n)
Return a structure used for the evaluation of sparse Hessian matrix in COO-format.
"""
struct HessSparseCOO{Tv, Ti} <: AbstractHess
  H::Symmetric{Tv, SparseMatrixCOO{Tv, Ti}}
end

function HessSparseCOO(nlp::AbstractNLPModel{T, S}, n) where {T, S}
  rows, cols = hess_structure(nlp)
  vals = S(undef, nlp.meta.nnzh)
  H = Symmetric(SparseMatrixCOO(n, n, rows, cols, vals), :L)
  return HessSparseCOO(H)
end

"""
    HessOp(::AbstractNLPModel{T,S}, n)
Return a structure used for the evaluation of the Hessian matrix as an operator.
"""
mutable struct HessOp{S, Op} <: AbstractHess
  Hv::S
  x_op::S
  H::Op
  function HessOp(nlp::AbstractNLPModel{T, S}, n) where {T, S}
    H = LinearOperator{T}(n, n, true, true, v -> v, v -> v, v -> v)
    x_op = copy(nlp.meta.x0)
    Hv = S(undef, n)
    H = hess_op!(nlp, x_op, Hv)
    return new{S, typeof(H)}(Hv, x_op, H)
  end
end

"""
    HessGaussNewtonOp(::AbstractNLSModel{T,S}, n)
Return a structure used for the evaluation of the Hessian matrix as an operator.
"""
mutable struct HessGaussNewtonOp{S, Op} <: AbstractHess
  Jv::S
  Jtv::S
  x_op::S
  H::Op
  function HessGaussNewtonOp(nls::AbstractNLSModel{T, S}, n) where {T, S}
    Jv, Jtv = S(undef, nls.nls_meta.nequ), S(undef, n)
    x_op = copy(nls.meta.x0)
    Jx = jac_op_residual!(nls, x_op, Jv, Jtv)
    H = Jx' * Jx
    return new{S, typeof(H)}(Jv, Jtv, x_op, H)
  end
end

export HessDense, HessSparse, HessSparseCOO, HessOp, HessGaussNewtonOp

"""
    init(::Type{Hess}, nlp::AbstractNLPModel{T,S}, n)

Return the hessian structure `Hess` and its composite type.
"""
function init(::Type{Hess}, nlp::AbstractNLPModel{T, S}, n) where {Hess <: AbstractHess, T, S}
  Hstruct = Hess(nlp, n)
  return Hstruct, typeof(Hstruct)
end

"""
    hessian!(workspace, nlp, x)

Return the Hessian matrix of `nlp` at `x` in-place with memory update of `workspace`.
"""
function hessian!(workspace::AbstractHess, nlp, x) end

function hessian!(workspace::HessDense, nlp, x)
  workspace.H .= Matrix(hess(nlp, x))
  return workspace.H
end

function hessian!(workspace::HessOp, nlp, x)
  workspace.x_op .= x
  return workspace.H
end

function hessian!(workspace::HessGaussNewtonOp, nlp, x)
  workspace.x_op .= x
  return workspace.H
end

function hessian!(workspace::HessSparse, nlp, x)
  hess_coord!(nlp, x, workspace.vals)
  n = nlp.meta.nvar
  workspace.H .= Symmetric(sparse(workspace.rows, workspace.cols, workspace.vals, n, n), :L)
  return workspace.H
end

function hessian!(workspace::HessSparseCOO, nlp, x)
  hess_coord!(nlp, x, workspace.H.data.vals)
  return workspace.H
end
