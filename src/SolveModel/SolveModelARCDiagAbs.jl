function solve_modelARCDiagAbs(H, g, nlp_stop, PData::PDataFact, α::T) where {T}
    #printsyled("on est dans solve_modelARCDiagAbs \n", color = :cyan)
    # Solve the ARC subproblem once diagonalized into Δ using the norm |Δ|

    # Setup the subproblem to recover identity regularizing matrix
    # Should this be in pre-process so that PData.λ is indeed used as updated?
    # No No! Use PData.success to initialize PData.λ when true, and use it otherwise
    #
    # TO TEST to confirm it avoids ascent directions
    #
    ḡ = PData.g̃
    #printsyled("on a ḡ = $ḡ\n", color = :yellow)
    n_g = norm(ḡ)
    #printsyled("on a n_g = $n_g\n", color = :yellow)
    # ϵ =  1.0e-10
    ϵ = sqrt(eps(T)) * 100.0
    #printsyled("on a ϵ = $ϵ  et PData.Δ = $(PData.Δ)\n", color = :yellow)
    Γ2 = max(maximum(abs.(PData.Δ)), ϵ)
    #printsyled("on a Γ2 = $Γ2 \n", color = :yellow)
    Γ = sqrt(Γ2)
    #printsyled("on a Γ = $Γ \n", color = :yellow)
    Δ̄ = Γ .\ PData.Δ ./ Γ
    #printsyled("on a Δ̄  = $Δ̄  \n", color = :yellow)
    Δ = PData.Δ
    #printsyled("on a Δ = $Δ \n", color = :yellow)
    λmin = 0.0 #  Test the following
    #printsyled("λmin = $λmin \n", color = :yellow)
    if PData.success # ensure
        #printsyled("PData.success = $(PData.success) \n", color = :yellow)
        l_m, = findmin(Δ̄)
        #printsyled("l_m = $l_m \n", color = :yellow)
        λ = max(-l_m, 0.0)
        #printsyled("λ = $λ \n", color = :yellow)
        λmin = max(ϵ, λ + ϵ * (1.0 + λ))
        #printsyled("λmin = $λmin \n", color = :yellow)
    else
        #printsyled("PData.success = $(PData.success) \n", color = :yellow)
        λmin = PData.λ
    end
    d̄ = -(Δ .+ λmin * Γ2) .\ ḡ
    #printsyled("d̄ = $d̄ \n", color = :yellow)
    seuil_bar = norm(d̄) / α

    #printsyled("on est a seuil_bar = $seuil_bar \n", color = :yellow)
    # Solve the subproblem (Δ̄ + λ I) d̄ = -ḡ such that λ = ||d̄||/α
    d̄, λ = solve_diag(λmin, Δ, ḡ, seuil_bar, α, ϵ, M = Γ2)
    # Transform back d_s into d
    d̃ = d̄
    #d̃ = Γ .\ d̄
    d = AInv(PData, d̃)

    #try assert((PData.g̃ + 0.5*PData.Δ .* d̃)⋅d̃ <= 0.0)  catch  @bp  end
    #println("*******SolveModelDiagAbs:  PData.g̃⋅d̃ = $(PData.g̃⋅d̃), 0.5 d̃'PData.Δd̃ = $(0.5*(PData.Δ .* d̃)⋅d̃)")

    PData.λ = λ

    return d, λ
end
