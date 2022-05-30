function solve_modelARCDiagAbs(H, g, nlp_stop, PData::PDataFact, α::T) where {T}
    # Solve the ARC subproblem once diagonalized into Δ using the norm |Δ|

    # Setup the subproblem to recover identity regularizing matrix
    # Should this be in pre-process so that PData.λ is indeed used as updated?
    # No No! Use PData.success to initialize PData.λ when true, and use it otherwise
    #
    # TO TEST to confirm it avoids ascent directions
    #
    ḡ = PData.g̃
    ϵ = sqrt(eps(T)) * 100.0
    Γ2 = max(maximum(abs.(PData.Δ)), ϵ)
    Γ = sqrt(Γ2)
    Δ̄ = Γ .\ PData.Δ ./ Γ
    Δ = PData.Δ
    λmin = 0.0 #  Test the following
    if PData.success
        l_m, = findmin(Δ̄)
        λ = max(-l_m, 0.0)
        λmin = max(ϵ, λ + ϵ * (1.0 + λ))
    else
        λmin = PData.λ
    end
    d̄ = -(Δ .+ λmin * Γ2) .\ ḡ
    seuil_bar = norm(d̄) / α

    # Solve the subproblem (Δ̄ + λ I) d̄ = -ḡ such that λ = ||d̄||/α
    d̄, λ = solve_diag(λmin, Δ, ḡ, seuil_bar, α, ϵ, M = Γ2)
    # Transform back d_s into d
    d̃ = d̄
    #d̃ = Γ .\ d̄
    d = AInv(PData, d̃)

    PData.λ = λ

    return d, λ
end
