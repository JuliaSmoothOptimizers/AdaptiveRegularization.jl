function solve_modelARCDiag(nlp_stop, PData :: PDataFact, α:: T) where T
    # Solve the ARC subproblem once diagonalized into Δ
    # printstyled("On est dans solve_modelARCDiag ⇊ \n", color = :red)
    # Setup the problem
    # @show eltype(α)
    # M = ones(PData.Δ)
    # M = fill(1, size(PData.Δ))
    M = ones(T, size(PData.Δ))
    # @show eltype(M)
    # printstyled("On a M ↢ \n", color = :red)

    # ϵ = 1.0e-10 * (1.0 + PData.λ)
    ϵ1 = sqrt(eps(T)) ./ T(100.0)
    # @show eltype(ϵ1)
    # printstyled("On a ϵ1 ↢  \n", color = :red)
    ϵ = T.(ϵ1 .* (1.0 .+ PData.λ))
    # @show eltype(ϵ)
    # printstyled("On a ϵ ↢  \n", color = :red)
    λ = max(ϵ, PData.λ .+ ϵ)
    # @show eltype(λ)
    # printstyled("On a λ ↢  \n", color = :red)

    d̃ = -(PData.Δ .+ λ .* M) .\ PData.g̃
    # @show eltype(d̃)
    # printstyled("On a d̃ ↢  \n", color = :red)
    normd̃ = sqrt(d̃⋅d̃)
    # @show eltype(normd̃)
    # printstyled("On a normd̃ ↢  \n", color = :red)
    seuil = normd̃ ./ α
    # @show eltype(seuil)
    # printstyled("On a seuil ↢  \n", color = :red)


    # Solve the subproblem (Δ + λ I) d̃ = -g̃ such that λ = ||d̃||/α
    d̃,λ = solve_diag(λ, PData.Δ, PData.g̃, seuil, α, ϵ)
    # @show eltype(d̃)
    # @show eltype(λ)
    # printstyled("On a d̃, λ ↢  \n", color = :red)

    PData.λ = λ
    # Transform back d̃ into d

    d = AInv(PData,d̃)
    # @show eltype(d)

    return d, NaN * rand(length(d)), λ
end
