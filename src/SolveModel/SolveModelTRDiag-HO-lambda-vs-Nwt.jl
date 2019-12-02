export solve_modelTRDiag_HO_vs_Nwt_λ

"""
If the Newton direction is accepted and the high order correction lies within
a bigger trust region then we use the high order correction.
"""
function solve_modelTRDiag_HO_vs_Nwt_λ(nlp_stop, PData :: PDataFact, δ:: T; ho_correction :: Symbol = :Shamanskii, fact = 2.0, λfact = 1.0, nwt_res_fact = 0.25) where T
    # Solve the TR subproblem once diagonalized into Δ using the norm |Δ|
    # Setup the problem
    # printstyled("On est dans solve_model_TRDiag_HO_vs_Nwt_λ  ↓ \n", color = :red)
    nlp_at_x = nlp_stop.current_state
    M = ones(T, size(PData.Δ))
    ϵ = sqrt(eps(T)) / T(100.0)
    ϵ2 = ϵ .* (T(1.0) .+ PData.λ)
    global dHO = nothing


    if PData.success # take care of eventual hard case and Newton's direction interior (λ = 0)
        # (PData.Δ + PData.λ * M) ⪰ 0
        # println("on a PData.succes = $(PData.success)")
        λ = max(ϵ2, ϵ2 + PData.λ) # to make sure (PData.Δ + λ * M) ≻ 0
        # λ = PData.λ

        d̃ = -(PData.Δ .+ λ * M) .\ PData.g̃ # Ajouter Shamanskii ici!
        normd̃ = sqrt(d̃⋅d̃)
        if normd̃ < δ
            if PData.λ == 0.0 # Newton's direction
                λ = PData.λ
                # println(" Newton's direction inside the region")
                #  d̃ is the Newton's direction, nothing more to do
            else              # hard case
                # println(" hard case")
                bidon, i = findmin(PData.Δ)
                d̃[i] = 0.0; d̃[i] = - sign(PData.g̃[i]) * sqrt(δ^2 - d̃⋅d̃)
            end
        else
            d̃,λ = solve_diagTR(λ, PData.Δ, PData.g̃, δ, ϵ)
        end
    else # hard case impossible, λ already > λ_min
        # println("on est dans le cas hard cases impossible")
        # Solve the subproblem (Δ + λ I) d̃ = -g̃ such that ||d̃|| = δ
        d̃, λ = solve_diagTR(PData.λ, PData.Δ, PData.g̃, δ, ϵ)
    end

    PData.λ = λ

    # Transform back d̃ into d

    d = AInv(PData, d̃)
    # @show λfact
    if λ <= λfact
        # printstyled("on a λ < λfact \n", color = :yellow)
        # @show ho_correction
        dHO = eval(ho_correction)(nlp_stop, PData, d, PData.g̃)
        # @show dHO
        # @show d
        # @show dHO == d
        if dHO == d
            # println("on est ici")
            return d, NaN * rand(length(d)), λ
        end
        # printstyled("on a dHO \n", color = :yellow)
        nwt_residual = (-(nlp_at_x.gx + 0.5 * nlp_at_x.Hx * d)⋅d)
        # printstyled("on a nwt_residual \n", color = :yellow)
        if (norm(dHO) < 2.0 .* δ) && ((-(nlp_at_x.gx + 0.5 * nlp_at_x.Hx * dHO)⋅dHO) >= nwt_res_fact .* nwt_residual)
            return dHO, dHO, λ
        else
            return d, dHO, λ
        end
    end

    #try assert((PData.g̃ + 0.5*PData.Δ .* d̃)⋅d̃ <= 0.0)  catch  @bp  end
    return d, NaN * rand(length(d)), λ
end
