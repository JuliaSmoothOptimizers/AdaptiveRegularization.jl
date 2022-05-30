export solve_modelTRDiagAbs
function solve_modelTRDiagAbs(H, g, nlp_stop, PData::PDataFact, δ::T) where {T}
    # Solve the TR subproblem once diagonalized into Δ using the norm |Δ|
    # Setup the subproblem to recover identity regularizing matrix
    # Should this be in pre-process so that PData.λ is indeed used as updated?

    ḡ = PData.g̃
    n_g = norm(ḡ)
    ϵ = sqrt(eps(T)) / T(100.0) # *n_g
    Γ2 = max(maximum(abs.(PData.Δ)), ϵ)
    Γ = sqrt(Γ2)

    Δ̄ = Γ .\ PData.Δ ./ Γ
    Δ = PData.Δ
    if PData.success    # initialize λmin, take care of hard case and Newton's direction interior
        l_m, i = findmin(Δ̄)
        λ = max(-l_m, 0.0)
        λmin = max(ϵ, λ .+ ϵ * (1.0 .+ λ)) #λ + ϵ #max(λ*(1+ϵ), ϵ)

        d̄ = -(Δ .+ λmin .* Γ2) .\ ḡ
        normd_s = sqrt(d̄ ⋅ d̄)
        if normd_s < δ
            if λ != 0.0  # λ == 0 is Newton's direction, nothing to do...
                d̄[i] = 0.0
                d̄[i] = -sign(PData.g̃[i]) * sqrt(δ^2 - d̄ ⋅ d̄)
            end
        else
            d̄, λ = solve_diagTR(λmin, Δ, ḡ, δ, ϵ, M = Γ2)
        end
    else # hard case impossible, λ already > λ_min
        # Solve the subproblem (Δ + λ I) d_s = -g_s such that δ = ||d_s||
        d̄, λ = solve_diagTR(PData.λ, Δ, ḡ, δ, ϵ, M = Γ2)
    end

    PData.λ = λ

    # Transform back d̄ into d
    d̃ = d̄ #Γ .\ d̄
    d = AInv(PData, d̃)
    return d, λ
end
