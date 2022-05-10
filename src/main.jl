export TRARC

struct TRARCWorkspace{T,S,Hess}
    xt::S
    xtnext::S
    d::S
    ∇f::S
    ∇fnext::S
    Hstruct::Hess
    function TRARCWorkspace(nlp::AbstractNLPModel{T, S}, ::Type{Hess}, n) where {T,S,Hess}
        return new{T,S,Hess}(
            S(undef, n), # xt
            S(undef, n), # xtnext
            S(undef, n), # d
            S(undef, n), # ∇f
            S(undef, n), # ∇fnext
            Hess(nlp, n),
        )
    end
end

function TRARC(
    nlp_stop::NLPStopping{Pb,M,SRC,NLPAtX{T,S},MStp,LoS};
    TR::TrustRegion = TrustRegion(T(10.0)),
    hess_type::Type{Hess} = HessDense,
    pdata_type::Type{ParamData} = PDataLDLt,
    solve_model::Function = solve_modelTRDiag,
    robust::Bool = true,
    verbose::Bool = false,
    kwargs...,
) where {Pb,M,SRC,MStp,LoS,S,T,Hess,ParamData}
    nlp, nlp_at_x = nlp_stop.pb, nlp_stop.current_state

    PData = ParamData(S, T, nlp.meta.nvar; kwargs...)
    workspace = TRARCWorkspace(nlp, Hess, nlp.meta.nvar)
    xt, xtnext, d, ∇f, ∇fnext =
        workspace.xt, workspace.xtnext, workspace.d, workspace.∇f, workspace.∇fnext

    α = TR.α₀  # initial Trust Region size
    max_unsuccinarow = TR.max_unsuccinarow
    acceptance_threshold = TR.acceptance_threshold
    increase_threshold = TR.increase_threshold
    reduce_threshold = TR.reduce_threshold

    xt .= nlp_at_x.x
    ft, ∇f = objgrad!(nlp, xt, ∇f)
    ftnext = ft
    norm_∇f = norm(∇f)
    nlp_stop.meta.optimality0 = norm_∇f

    OK = update_and_start!(nlp_stop, x = xt, fx = ft, gx = ∇f)
    !OK && Stopping.update!(
        nlp_at_x,
        Hx = hessian!(workspace.Hstruct, nlp, xt),
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
        max_hprod = nlp_stop.meta.max_cntrs[:neval_hprod]
        PData = preprocess(PData, nlp_at_x.Hx, ∇f, neval_hprod(nlp), max_hprod, α)

        if ~PData.OK
            @warn("Something wrong with PData")
            return nlp_stop
        end

        success = false
        while !success & (unsuccinarow < max_unsuccinarow)
            d, λ = solve_model(nlp_at_x.Hx, ∇f, nlp_stop, PData, α) # Est-ce que d et λ ne sont pas dans PData ?

            Δq = -(∇f + 0.5 * nlp_at_x.Hx * d) ⋅ d

            if Δq < 0.0
                println("*******   Ascent direction in SolveModel: Δq = $Δq")
                println(
                    "  g⋅d = $(∇f⋅d), 0.5 d'Hd = $(0.5*(nlp_at_x.Hx*d)⋅d)  α = $α  λ = $λ",
                )
                # cond issue with H?
                return nlp_stop
            end
            slope = ∇f ⋅ d
            xtnext .= xt .+ d
            ftnext = obj(nlp, xtnext)
            Δf = ft - ftnext

            iter += 1

            r, good_grad, ∇fnext =
                compute_r(nlp, ft, Δf, Δq, slope, d, xtnext, ∇fnext, robust)

            if r < acceptance_threshold # unsucessful
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
                    grad!(nlp, xt, ∇f)
                end

                verysucc += 1
                if r > increase_threshold # very sucessful
                    α = increase(PData, α, TR)
                    verbose &&
                        @info log_row(Any[iter, ft, stop_norm(∇f), λ, "V", α, norm(d), Δq])
                else # sucessful
                    if r < reduce_threshold
                        α = decrease(PData, α, TR)
                    end
                    verbose &&
                        @info log_row(Any[iter, ft, stop_norm(∇f), λ, "S", α, norm(d), Δq])
                    succ += 1
                end
            end
        end # while !success

        nlp_stop.meta.nb_of_stop = iter
        OK = update_and_stop!(nlp_stop, x = xt, fx = ft, gx = ∇f)
        success && Stopping.update!(nlp_at_x, Hx = hessian!(workspace.Hstruct, nlp, xt))
    end # while !OK

    return nlp_stop
end
