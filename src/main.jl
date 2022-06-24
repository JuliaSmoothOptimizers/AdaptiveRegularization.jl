"""
    TRARCWorkspace(nlp, ::Type{Hess}, n)
Pre-allocate the memory used during the [`TRARC`](@ref) call for the problem `nlp` of size `n`.
The possible values for `Hess` are: `HessDense`, `HessSparse`, `HessSparseCOO`, `HessOp`.
Return a `TRARCWorkspace` structure.
"""
struct TRARCWorkspace{T,S,Hess}
    xt::S
    xtnext::S
    d::S
    ∇f::S
    ∇fnext::S
    Hstruct::Hess
    Fx::S
    function TRARCWorkspace(nlp::AbstractNLPModel{T,S}, ::Type{Hess}) where {T,S,Hess}
        n = nlp.meta.nvar
        return new{T,S,Hess}(
            S(undef, n), # xt
            S(undef, n), # xtnext
            S(undef, n), # d
            S(undef, n), # ∇f
            S(undef, n), # ∇fnext
            Hess(nlp, n),
            S(undef, 0),
        )
    end
    function TRARCWorkspace(nlp::AbstractNLSModel{T,S}, ::Type{Hess}) where {T,S,Hess}
        n = nlp.meta.nvar
        return new{T,S,Hess}(
            S(undef, n), # xt
            S(undef, n), # xtnext
            S(undef, n), # d
            S(undef, n), # ∇f
            S(undef, n), # ∇fnext
            Hess(nlp, n),
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

function preprocess(stp::NLPStopping, PData::TPData, workspace::TRARCWorkspace, ∇f, norm_∇f, α)
    max_hprod = stp.meta.max_cntrs[:neval_hprod]
    Hx = stp.current_state.Hx
    PData = preprocess(PData, Hx, ∇f, norm_∇f, neval_hprod(stp.pb), max_hprod, α)
    return PData
end

function compute_direction(stp::NLPStopping, PData::TPData, workspace::TRARCWorkspace, ∇f, norm_∇f, α, solve_model)
    max_hprod = stp.meta.max_cntrs[:neval_hprod]
    Hx = stp.current_state.Hx
    return solve_model(PData, Hx, ∇f, norm_∇f, neval_hprod(stp.pb), max_hprod, α)
end

function compute_direction(stp::NLPStopping, PData::PDataIterLS, workspace::TRARCWorkspace, ∇f, norm_∇f, α, solve_model)
    max_prod = stp.meta.max_cntrs[:neval_jprod_residual]
    Jx = jac_op_residual(stp.pb, workspace.xt)
    Fx = workspace.Fx
    return solve_model(PData, Jx, Fx, norm_∇f, neval_jprod_residual(stp.pb), max_prod, α)
end

function compute_direction(stp::NLPStopping, PData::PDataIterLS, workspace::TRARCWorkspace{T,S,Hess}, ∇f, norm_∇f, α, solve_model) where {T,S,Hess <: HessGaussNewtonOp}
    max_prod = stp.meta.max_cntrs[:neval_jprod_residual]
    Jx = jac_op_residual!(stp.pb, workspace.xt, workspace.Hstruct.Jv, workspace.Hstruct.Jtv)
    Fx = workspace.Fx
    return solve_model(PData, Jx, Fx, norm_∇f, neval_jprod_residual(stp.pb), max_prod, α)
end

function hessian!(workspace::TRARCWorkspace, nlp, x)
    return hessian!(workspace.Hstruct, nlp, x)
end

function TRARC(
    nlp_stop::NLPStopping{Pb,M,SRC,NLPAtX{T,S},MStp,LoS};
    TR::TrustRegion = TrustRegion(T(10.0)),
    hess_type::Type{Hess} = HessOp,
    pdata_type::Type{ParamData} = PDataKARC,
    kwargs...,
) where {Pb,M,SRC,MStp,LoS,S,T,Hess,ParamData}
    nlp = nlp_stop.pb

    if ParamData == PDataNLSST
        PData = PDataNLSST(S, T, nlp.meta.nvar, nlp.nls_meta.nequ; kwargs...)
    else
        PData = ParamData(S, T, nlp.meta.nvar; kwargs...)
    end
    workspace = TRARCWorkspace(nlp, Hess)
    return TRARC(nlp_stop, PData, workspace, TR; kwargs...)
end

function TRARC(
    nlp_stop::NLPStopping{Pb,M,SRC,NLPAtX{T,S},MStp,LoS},
    PData::ParamData,
    workspace::TRARCWorkspace{T,S,Hess},
    TR::TrustRegion;
    solve_model::Function = solve_modelKARC,
    robust::Bool = true,
    verbose::Bool = false,
    kwargs...,
) where {Pb,M,SRC,MStp,LoS,S,T,Hess,ParamData}
    nlp, nlp_at_x = nlp_stop.pb, nlp_stop.current_state
    xt, xtnext, d, ∇f, ∇fnext =
        workspace.xt, workspace.xtnext, workspace.d, workspace.∇f, workspace.∇fnext

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

    OK = update_and_start!(nlp_stop, x = xt, fx = ft, gx = ∇f)
    !OK && Stopping.update!(
        nlp_at_x,
        Hx = hessian!(workspace, nlp, xt),
        convert = true,
    )

    iter = 0 # counter different than stop count
    succ, unsucc, verysucc, unsuccinarow = 0, 0, 0, 0

    verbose && @info log_header(
        [:iter, :f, :nrm_g, :λ, :status, :α, :nrm_dtr, :f_p_dTR, :ΔqN],
        [Int64, T, T, T, String, T, T, T],
    )
    verbose && @info log_row(Any[iter, ft, norm_∇f, 0.0, "First iteration", α])

    while !OK
        PData = preprocess(nlp_stop, PData, workspace, ∇f, norm_∇f, α)

        if ~PData.OK
            @warn("Something wrong with PData")
            return nlp_stop
        end

        success = false
        while !success & (unsuccinarow < max_unsuccinarow)
            d, λ = compute_direction(nlp_stop, PData, workspace, ∇f, norm_∇f, α, solve_model)

            slope = ∇f ⋅ d
            Δq = -(∇f + 0.5 * (nlp_at_x.Hx * d)) ⋅ d

            xtnext .= xt .+ d
            ftnext = obj(nlp, xtnext, workspace)
            Δf = ft - ftnext

            iter += 1

            r, good_grad, ∇fnext =
                compute_r(nlp, ft, Δf, Δq, slope, d, xtnext, workspace, robust)

            if Δq < 0.0 # very unsucessful
                verbose && @info log_row(Any[iter, ft, norm_∇f, λ, "VU", α, norm(d), Δq])
                unsucc += 1
                unsuccinarow += 1
                η = (1 - acceptance_threshold) / 10 # ∈ (acceptance_threshold, 1)
                qksk = ft + slope + ((nlp_at_x.Hx * d) ⋅ d) / 2
                αbad = (1 - η) * slope / ((1 - η) * (ft + slope) + η * qksk - ftnext)
                α = min(decrease(PData, α, TR), max(TR.large_decrease_factor, αbad) * α)
            elseif r < acceptance_threshold # unsucessful
                verbose && @info log_row(Any[iter, ft, norm_∇f, λ, "U", α, norm(d), Δq])
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
                    verbose && @info log_row(Any[iter, ft, norm_∇f, λ, "V", α, norm(d), Δq])
                else # sucessful
                    if r < reduce_threshold
                        α = decrease(PData, α, TR)
                    end
                    verbose && @info log_row(Any[iter, ft, norm_∇f, λ, "S", α, norm(d), Δq])
                    succ += 1
                end
            end
        end # while !success

        nlp_stop.meta.nb_of_stop = iter
        OK = update_and_stop!(nlp_stop, x = xt, fx = ft, gx = ∇f)
        success && Stopping.update!(nlp_at_x, Hx = hessian!(workspace, nlp, xt))
    end # while !OK

    return nlp_stop
end
