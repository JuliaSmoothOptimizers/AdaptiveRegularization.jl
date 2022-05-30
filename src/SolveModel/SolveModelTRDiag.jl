export solve_modelTRDiag

function solve_modelTRDiag(H, g, nlp_stop, PData::PDataFact, δ::T) where {T}
    # Solve the TR subproblem once diagonalized into Δ using the norm |Δ|
    M = ones(T, size(PData.Δ))
    ϵ = sqrt(eps(T)) / T(100.0)
    ϵ2 = ϵ .* (T(1.0) .+ PData.λ)
    global dHO = nothing


    if PData.success # take care of eventual hard case and Newton's direction interior (λ = 0)
        λ = max(ϵ2, ϵ2 + PData.λ) # to make sure (PData.Δ + λ * M) ≻ 0

        d̃ = -(PData.Δ .+ λ * M) .\ PData.g̃
        normd̃ = sqrt(d̃ ⋅ d̃)
        if normd̃ < δ
            if PData.λ == 0.0 # Newton's direction
                λ = PData.λ
                #  d̃ is the Newton's direction, nothing more to do
            else              # hard case
                _, i = findmin(PData.Δ)
                d̃[i] = 0.0
                d̃[i] = -sign(PData.g̃[i]) * sqrt(δ^2 - d̃ ⋅ d̃)
            end
        else
            d̃, λ = solve_diagTR(λ, PData.Δ, PData.g̃, δ, ϵ)
        end
    else # hard case impossible, λ already > λ_min
        # Solve the subproblem (Δ + λ I) d̃ = -g̃ such that ||d̃|| = δ
        d̃, λ = solve_diagTR(PData.λ, PData.Δ, PData.g̃, δ, ϵ)
    end

    PData.λ = λ

    # Transform back d̃ into d
    d = AInv(PData, d̃)

    return d, λ
end
