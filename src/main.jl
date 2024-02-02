function preprocess!(
  stp::NLPStopping{Pb, M, SRC, NLPAtX{Score, T, S}, MStp, LoS},
  PData::TPData{T},
  workspace::TRARCWorkspace{T, S, Hess},
  ∇f::S,
  norm_∇f::T,
  α::T,
) where {Pb, M, SRC, MStp, LoS, Score, S, T, Hess}
  max_hprod = stp.meta.max_cntrs[:neval_hprod]
  Hx = get_hess(workspace.Hstruct)
  preprocess!(PData, Hx, ∇f, norm_∇f, neval_hprod(stp.pb), max_hprod, α)
  return PData
end

function compute_direction(
  stp::NLPStopping{Pb, M, SRC, NLPAtX{Score, T, S}, MStp, LoS},
  PData::TPData{T},
  workspace::TRARCWorkspace{T, S, Hess},
  ∇f::S,
  norm_∇f::T,
  α::T,
) where {Pb, M, SRC, MStp, LoS, Score, S, T, Hess}
  max_hprod = stp.meta.max_cntrs[:neval_hprod]
  Hx = get_hess(workspace.Hstruct)
  solve_model!(PData, Hx, ∇f, norm_∇f, neval_hprod(stp.pb), max_hprod, α)
  return PData.d, PData.λ
end

function compute_direction(
  stp::NLPStopping{Pb, M, SRC, NLPAtX{Score, T, S}, MStp, LoS},
  PData::PDataIterLS{T},
  workspace::TRARCWorkspace{T, S, Hess},
  ∇f::S,
  norm_∇f::T,
  α::T,
) where {Pb, M, SRC, MStp, LoS, Score, S, T, Hess}
  max_prod = stp.meta.max_cntrs[:neval_jprod_residual]
  Jx = jac_op_residual(stp.pb, workspace.xt)
  Fx = workspace.Fx
  solve_model!(PData, Jx, Fx, norm_∇f, neval_jprod_residual(stp.pb), max_prod, α)
  return PData.d, PData.λ
end

function compute_direction(
  stp::NLPStopping{Pb, M, SRC, NLPAtX{Score, T, S}, MStp, LoS},
  PData::PDataIterLS{T},
  workspace::TRARCWorkspace{T, S, Hess},
  ∇f::S,
  norm_∇f::T,
  α::T,
) where {Pb, M, SRC, MStp, LoS, Score, S, T, Hess <: HessGaussNewtonOp}
  max_prod = stp.meta.max_cntrs[:neval_jprod_residual]
  Jx = jac_op_residual!(stp.pb, workspace.xt, workspace.Hstruct.Jv, workspace.Hstruct.Jtv)
  Fx = workspace.Fx
  solve_model!(PData, Jx, Fx, norm_∇f, neval_jprod_residual(stp.pb), max_prod, α)
  return PData.d, PData.λ
end

function hessian!(workspace::TRARCWorkspace, nlp, x)
  return hessian!(workspace.Hstruct, nlp, x)
end

"""
    compute_Δq(Hx, d, ∇f)

Update `Δq = -(∇f + 0.5 * (Hx * d)) ⋅ d` in-place.
"""
function compute_Δq(workspace, Hx, d, ∇f)
  mul!(workspace.Hd, Hx, d)
  workspace.Hd .*= 1 // 2
  workspace.Hd .+= ∇f
  return -dot(workspace.Hd, d)
end

function SolverCore.solve!(
  solver::TRARCSolver{T, V},
  nlp::AbstractNLPModel{T, V},
  stats::GenericExecutionStats{T, V};
  x::V = nlp.meta.x0,
  atol::T = √eps(T),
  rtol::T = √eps(T),
  max_time::Float64 = 300.0,
  kwargs...,
) where {T, V}
  stp = solver.stp
  stp.meta.atol = atol
  stp.meta.rtol = rtol
  stp.meta.max_time = max_time
  if x != stp.current_state.x
    set_x!(stp.current_state, x)
    grad!(nlp, x, solver.workspace)
    set_gx!(stp.current_state, solver.workspace.∇f)
    set_res!(stp.current_state, stp.current_state.gx)
    # we would also need to reinit the `tol_check` function
  end
  return SolverCore.solve!(solver, stp, stats; kwargs...)
end

# Main algorithm
function SolverCore.solve!(
  solver::TRARCSolver{T, S},
  nlp_stop::NLPStopping{Pb, M, SRC, NLPAtX{Score, T, S}, MStp, LoS},
  stats::GenericExecutionStats{T, S};
  robust::Bool = true,
  verbose::Integer = false,
  callback = (args...) -> nothing,
  kwargs...,
) where {Pb, M, SRC, MStp, LoS, Score, S, T}
  PData = solver.meta
  workspace = solver.workspace
  TR = solver.TR
  nlp, nlp_at_x = nlp_stop.pb, nlp_stop.current_state
  xt, xtnext, ∇f, ∇fnext = workspace.xt, workspace.xtnext, workspace.∇f, workspace.∇fnext
  d = workspace.d
  Hx = get_hess(workspace.Hstruct)
  reset!(stats)

  α = TR.α₀
  max_unsuccinarow = TR.max_unsuccinarow
  acceptance_threshold = TR.acceptance_threshold
  increase_threshold = TR.increase_threshold
  reduce_threshold = TR.reduce_threshold

  xt .= nlp_at_x.x
  ft, ∇f = objgrad!(nlp, xt, workspace)
  ftnext = ft
  norm_∇f = norm(∇f)
  nlp_stop.meta.optimality0 = norm_∇f

  set_x!(nlp_at_x, xt)
  set_fx!(nlp_at_x, ft)
  set_gx!(nlp_at_x, ∇f)
  OK = start!(nlp_stop)
  Hx = hessian!(workspace, nlp, xt)

  iter = 0 # counter different than stop count
  succ, unsucc, verysucc, unsuccinarow = 0, 0, 0, 0

  set_status!(stats, status_stopping_to_stats(nlp_stop))
  set_solution!(stats, nlp_at_x.x)
  set_objective!(stats, nlp_at_x.fx)
  set_dual_residual!(stats, nlp_at_x.current_score)
  set_iter!(stats, 0)
  set_time!(stats, nlp_at_x.current_time - nlp_stop.meta.start_time)

  verbose > 0 && @info log_header(
    [:iter, :f, :nrm_g, :λ, :status, :α, :nrm_dtr, :f_p_dTR, :ΔqN],
    [Int64, T, T, T, String, T, T, T],
  )
  verbose > 0 && @info log_row(Any[iter, ft, norm_∇f, 0.0, "First iteration", α])

  callback(nlp, solver, stats)

  while !OK && (stats.status != :user)
    preprocess!(nlp_stop, PData, workspace, ∇f, norm_∇f, α)

    if ~PData.OK
      @warn("Something wrong with PData")
      return nlp_stop
    end

    success = false
    while !success & (unsuccinarow < max_unsuccinarow)
      d, λ = compute_direction(nlp_stop, PData, workspace, ∇f, norm_∇f, α)

      slope = ∇f ⋅ d
      Δq = compute_Δq(workspace, Hx, d, ∇f)

      xtnext .= xt .+ d
      ftnext = obj(nlp, xtnext, workspace)
      Δf = ft - ftnext

      iter += 1

      r, good_grad, ∇fnext = compute_r(nlp, ft, Δf, Δq, slope, d, xtnext, workspace, robust)

      if Δq < 0.0 # very unsucessful
        verbose > 0 &&
          mod(iter, verbose) == 0 &&
          @info log_row(Any[iter, ft, norm_∇f, λ, "VU", α, norm(d), Δq])
        unsucc += 1
        unsuccinarow += 1
        η = (1 - acceptance_threshold) / 10 # ∈ (acceptance_threshold, 1)
        mul!(workspace.Hd, Hx, d)
        qksk = ft + slope + (workspace.Hd ⋅ d) / 2
        αbad = (1 - η) * slope / ((1 - η) * (ft + slope) + η * qksk - ftnext)
        α = min(decrease(PData, α, TR), max(TR.large_decrease_factor, αbad) * α)
      elseif r < acceptance_threshold # unsucessful
        verbose > 0 &&
          mod(iter, verbose) == 0 &&
          @info log_row(Any[iter, ft, norm_∇f, λ, "U", α, norm(d), Δq])
        unsucc += 1
        unsuccinarow += 1
        α = decrease(PData, α, TR)
      else
        success = true
        unsuccinarow = 0

        xt .= xtnext
        ft = ftnext
        if good_grad
          ∇f .= ∇fnext
        else
          ∇f = grad!(nlp, xt, workspace)
        end
        norm_∇f = norm(∇f)

        verysucc += 1
        if r > increase_threshold # very sucessful
          α = increase(PData, α, TR)
          verbose > 0 &&
            mod(iter, verbose) == 0 &&
            @info log_row(Any[iter, ft, norm_∇f, λ, "V", α, norm(d), Δq])
        else # sucessful
          if r < reduce_threshold
            α = decrease(PData, α, TR)
          end
          verbose > 0 &&
            mod(iter, verbose) == 0 &&
            @info log_row(Any[iter, ft, norm_∇f, λ, "S", α, norm(d), Δq])
          succ += 1
        end
      end
    end # while !success

    nlp_stop.meta.nb_of_stop = iter
    set_x!(nlp_at_x, xt)
    set_fx!(nlp_at_x, ft)
    set_gx!(nlp_at_x, ∇f)
    OK = stop!(nlp_stop)
    Hx = hessian!(workspace, nlp, xt)

    set_status!(stats, status_stopping_to_stats(nlp_stop))
    set_solution!(stats, nlp_at_x.x)
    set_objective!(stats, nlp_at_x.fx)
    set_dual_residual!(stats, nlp_at_x.current_score)
    set_iter!(stats, nlp_stop.meta.nb_of_stop)
    set_time!(stats, nlp_at_x.current_time - nlp_stop.meta.start_time)
    callback(nlp, solver, stats)
  end # while !OK

  stats
end
