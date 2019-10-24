export solve_modelTRDiag_HO_vs_Cauchy

"""
If the Newton direction is accepted and the high order correction lies within
a bigger trust region then we use the high order correction.
"""
function solve_modelTRDiag_HO_vs_Cauchy(nlp_stop, PData :: PDataFact, δ:: T; ho_correction :: Symbol = :Shamanskii) where T
    # Solve the TR subproblem once diagonalized into Δ using the norm |Δ|
    # Setup the problem
    # printstyled("On est dans solve_modelTRDiag_HO_vs_cauchy \n", color = :red)
    nlp_at_x = nlp_stop.current_state
    M = fill(T.(1.0), size(PData.Δ))
    ϵ = sqrt(eps(T)) / T(100.0)
    ϵ2 = T.(ϵ * (T(1.0) + PData.λ))
    global dHO = nothing
    global τ = nothing

    # log_header([:λ, :norm_dn, :norm_dho], [T, T, T, T, T, T])

    if PData.success # take care of eventual hard case and Newton's direction interior (λ = 0)
        # printstyled("on a PData.succes = $(PData.success) \n", color = :red)
        # (PData.Δ + PData.λ * M) ⪰ 0
        λ = max(ϵ2, PData.λ + ϵ2) # to make sure (PData.Δ + λ * M) ≻ 0

        d̃ = -(PData.Δ .+ λ * M) .\ PData.g̃
        normd̃ = sqrt(d̃⋅d̃)
        normg̃ = sqrt(PData.g̃⋅PData.g̃)
        if normd̃ < δ
            if PData.λ == 0.0 # Newton's direction
                # println(" Newton's direction inside the region")
                λ = PData.λ
                dN = AInv(PData, d̃)
                # @show dN
                dHO = eval(ho_correction)(nlp_stop, PData, dN, PData.g̃)
                # @show dHO
                gBg = (nlp_at_x.Hx * nlp_at_x.gx)' * nlp_at_x.gx
                ng = norm(nlp_at_x.gx)
                if gBg <= 0.0
                    τ = 1.0
                else
                    ng3 = ng^3
                    denom = δ .* gBg
                    τ = min(ng3 ./ denom, 1)
                end
                # @show τ
                dCauchy = -τ .* ((δ ./ ng) * nlp_at_x.gx)
                # @show dCauchy
                cauchy_decrease = -(nlp_at_x.gx + 0.5 * nlp_at_x.Hx * dCauchy)⋅dCauchy
                # @show cauchy_decrease
                if (norm(dHO) < 2.0 .* δ) && ((-(nlp_at_x.gx + 0.5 * nlp_at_x.Hx * dHO)⋅dHO) > cauchy_decrease) # && (norm(grad(nlp_stop.pb, nlp_at_x.x + dHO)) < norm(nlp_at_x.gx))
                    # printstyled("on prend dHO 🐣\n", color = :green)
                    return dHO, dHO, λ
                else
                    return dN, dHO, λ
                end
                # return dN, dHO, xdemi, λ
                # println(" Newton's direction inside the region")
                #  d̃ is the Newton's direction, nothing more to do
            else              # hard case
                # println(" hard case")
                bidon, i = findmin(PData.Δ)
                d̃[i] = T(0.0); d̃[i] = - sign(PData.g̃[i]) * sqrt(δ^2 - d̃⋅d̃)
            end
        else
            # println("success pas hard case")
            d̃,λ = solve_diagTR(λ, PData.Δ, PData.g̃, δ, ϵ)
        end
    else # hard case impossible, λ already > λ_min
        # println("on est dans le cas !(PData.success)")
        # Solve the subproblem (Δ + λ I) d̃ = -g̃ such that ||d̃|| = δ
        # printstyled("hard case impossible \n", color = :yellow)
        d̃, λ = solve_diagTR(PData.λ, PData.Δ, PData.g̃, δ, ϵ)
    end

    PData.λ = λ

    # Transform back d̃ into d
    d = AInv(PData, d̃)
    # printstyled("On a d après AInv \n", color = :red)
    return d, NaN * rand(length(d)), λ
    #try assert((PData.g̃ + 0.5*PData.Δ .* d̃)⋅d̃ <= 0.0)  catch  @bp  end
    # return d, λ
end
