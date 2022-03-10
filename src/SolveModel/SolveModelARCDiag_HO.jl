export solve_modelARCDiag_HO

function solve_modelARCDiag_HO(nlp_stop, PData :: PDataFact, α:: T; ho_correction :: Symbol = :Shamanskii, λfact = 1.0) where T
    # Solve the ARC subproblem once diagonalized into Δ
    # printstyled("On est dans solve_modelARCDiag ⇊ \n", color = :red)
    nlp_at_x = nlp_stop.current_state
    # Setup the problem
    # M = ones(PData.Δ)
    # M = fill(1, size(PData.Δ))
    M = ones(T, size(PData.Δ))
    # printstyled("On a M ↢ \n", color = :red)

    # ϵ = 1.0e-10 * (1.0 + PData.λ)
    ϵ1 = sqrt(eps(T)) / 100.0
    # printstyled("On a ϵ1 ↢  \n", color = :red)
    ϵ = T.(ϵ1 * (1.0 + PData.λ))
    # printstyled("On a ϵ ↢  \n", color = :red)
    λ = max(ϵ,PData.λ+ϵ)
    # printstyled("On a λ ↢  \n", color = :red)

    d̃ = -(PData.Δ+λ*M) .\ PData.g̃
    # printstyled("On a d̃ ↢  \n", color = :red)
    normd̃ = sqrt(d̃⋅d̃)
    # printstyled("On a normd̃ ↢  \n", color = :red)
    seuil = normd̃/α
    # printstyled("On a seuil ↢  \n", color = :red)


    # Solve the subproblem (Δ + λ I) d̃ = -g̃ such that λ = ||d̃||/α
    d̃,λ = solve_diag(λ,PData.Δ,PData.g̃,seuil,α,ϵ)
    # printstyled("On a d̃, λ ↢  \n", color = :red)

    PData.λ = λ
    # Transform back d̃ into d

    d = AInv(PData,d̃)

    if PData.λ <= λfact
        # println("ici ! 🍆")
        # @show λfact
        # @show ho_correction
        dHO = eval(ho_correction)(nlp_stop, PData, d, PData.g̃)
        # @show dHO
        # @show (norm(dHO) < 2.0 .* α)
        # @show ((-(nlp_at_x.gx + 0.5 .* nlp_at_x.Hx * dHO)⋅dHO) > 0.0)
        # @show nlp_at_x.gx
        # @show nlp_at_x.Hx
        # @show -(nlp_at_x.gx + 0.5 .* nlp_at_x.Hx * dHO)
        # @show (-(nlp_at_x.gx + 0.5 .* nlp_at_x.Hx * dHO)⋅dHO)

        if (norm(dHO) < 2.0 .* α) && ((-(nlp_at_x.gx + 0.5 * nlp_at_x.Hx * dHO)⋅dHO) > 0.0)
            # printstyled("on prend dHO 🐣\n", color = :green)
            return dHO, PData.λ
        else
            # printstyled("on prend d 🐲 \n", color = :green)
            return d, PData.λ
        end
    end

    return d, λ
end
