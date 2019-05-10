function solve_modelARCDiag(PData :: PDataFact, α:: T) where T
    # Solve the ARC subproblem once diagonalized into Δ

    # Setup the problem
    # M = ones(PData.Δ)
    M = fill(1, size(PData.Δ))

    # ϵ = 1.0e-10 * (1.0 + PData.λ)
    ϵ1 = sqrt(eps(T)) / 100.0
    ϵ = T.(ϵ1 * (1.0 + PData.λ))
    λ = max(ϵ,PData.λ+ϵ)

    d̃ = -(PData.Δ+λ*M) .\ PData.g̃
    normd̃ = sqrt(d̃⋅d̃)
    seuil = normd̃/α


    # Solve the subproblem (Δ + λ I) d̃ = -g̃ such that λ = ||d̃||/α
    d̃,λ = solve_diag(λ,PData.Δ,PData.g̃,seuil,α,ϵ)

    PData.λ = λ
    # Transform back d̃ into d

    d = AInv(PData,d̃)

    return d, λ
end
