"""
    TRARCWorkspace(nlp, ::Type{Hess}, n)
Pre-allocate the memory used during the [`TRARC`](@ref) call for the problem `nlp` of size `n`.
The possible values for `Hess` are: `HessDense`, `HessSparse`, `HessSparseCOO`, `HessOp`.
Return a `TRARCWorkspace` structure.
"""
struct TRARCWorkspace{T, S, Hess}
  xt::S
  xtnext::S
  d::S
  ∇f::S
  ∇fnext::S
  Hstruct::Hess
  Hd::S
  Fx::S
  function TRARCWorkspace(nlp::AbstractNLPModel{T, S}, ::Type{Hess}) where {T, S, Hess}
    n = nlp.meta.nvar
    Hstruct, Htype = init(Hess, nlp, n)
    return new{T, S, Htype}(
      S(undef, n), # xt
      S(undef, n), # xtnext
      S(undef, n), # d
      S(undef, n), # ∇f
      S(undef, n), # ∇fnext
      Hstruct,
      S(undef, n), # H * d
      S(undef, 0),
    )
  end
  function TRARCWorkspace(nlp::AbstractNLSModel{T, S}, ::Type{Hess}) where {T, S, Hess}
    n = nlp.meta.nvar
    Hstruct, Htype = init(Hess, nlp, n)
    return new{T, S, Htype}(
      S(undef, n), # xt
      S(undef, n), # xtnext
      S(undef, n), # d
      S(undef, n), # ∇f
      S(undef, n), # ∇fnext
      Hstruct,
      S(undef, n), # H * d
      S(undef, nlp.nls_meta.nequ),
    )
  end
end

# Redefined NLP Model API to use workspace
function NLPModels.objgrad!(nlp::AbstractNLPModel, x, workspace::TRARCWorkspace)
  return objgrad!(nlp, x, workspace.∇f)
end

function NLPModels.obj(nlp::AbstractNLPModel, x, workspace::TRARCWorkspace)
  return obj(nlp, x)
end

function NLPModels.grad!(nlp::AbstractNLPModel, x, workspace::TRARCWorkspace)
  return grad!(nlp, x, workspace.∇f)
end

function NLPModels.objgrad!(nls::AbstractNLSModel, x, workspace::TRARCWorkspace)
  increment!(nls, :neval_obj)
  increment!(nls, :neval_grad)
  Fx = residual!(nls, x, workspace.Fx)
  ∇f = jtprod_residual!(nls, x, Fx, workspace.∇f)
  return dot(Fx, Fx) / 2, ∇f
end

function NLPModels.obj(nls::AbstractNLSModel, x, workspace::TRARCWorkspace)
  increment!(nls, :neval_obj)
  Fx = residual!(nls, x, workspace.Fx)
  return dot(Fx, Fx) / 2
end

function NLPModels.grad!(nls::AbstractNLSModel, x, workspace::TRARCWorkspace)
  increment!(nls, :neval_grad)
  Fx = workspace.Fx
  return jtprod_residual!(nls, x, Fx, workspace.∇f)
end

"""
    TRARCSolver(nlp::AbstractNLPModel [, x0 = nlp.meta.x0]; kwargs...)
    TRARCSolver(stp::NLPStopping; kwargs...)

Structure regrouping all the structure used during the `TRARC` call. It returns a `TRARCSolver` structure.

# Arguments
The keyword arguments may include:

- `stp::NLPStopping`: `Stopping` structure for this algorithm workflow;
- `meta::ParamData`: see [`ParamData`](@ref);
- `workspace::TRARCWorkspace`: allocated space for the solver itself;
- `TR::TrustRegion`: trust-region parameters.

"""
mutable struct TRARCSolver{
  T,
  S,
  # Hess,
  Pb <: AbstractNLPModel{T, S},
} <: AbstractOptimizationSolver
  stp::NLPStopping{Pb}
  meta
  workspace # ::TRARCWorkspace{T, S, Hess}
  TR # ::TrustRegion
end

export TRARCSolver

function TRARCSolver(nlp::AbstractNLPModel; kwargs...)
  stp = NLPStopping(nlp; optimality_check = (pb, state) -> norm(state.gx), kwargs...)
  return TRARCSolver(stp; kwargs...)
end

function TRARCSolver(
  stp::NLPStopping{Pb, M, SRC, NLPAtX{Score, T, S}, MStp, LoS};
  TR::TrustRegion = TrustRegion(T(10.0)),
  hess_type::Type{Hess} = HessOp,
  pdata_type::Type{ParamData} = PDataKARC,
  kwargs...,
) where {Pb, M, SRC, MStp, LoS, Score, S, T, Hess, ParamData}
  nlp = stp.pb
  meta = if ParamData == PDataNLSST
    PDataNLSST(S, T, nlp.meta.nvar, nlp.nls_meta.nequ; kwargs...)
  else
    ParamData(S, T, nlp.meta.nvar; kwargs...)
  end
  workspace = TRARCWorkspace(nlp, Hess)
  return TRARCSolver(stp, meta, workspace, TR)
end

function SolverCore.reset!(solver::TRARCSolver)
  reinit!(solver.stp)
  solver
end

function SolverCore.reset!(solver::TRARCSolver, nlp::AbstractNLPModel)
  @assert nlp.meta.nvar == solver.stp.pb.meta.nvar
  @assert nlp.meta.ncon == solver.stp.pb.meta.ncon
  reset!(solver)
  solver.stp.pb = nlp
  reinit!(solver.stp)
  solver
end
