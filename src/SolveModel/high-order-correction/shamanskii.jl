export Shamanskii

function Shamanskii(
    nlp_stop,
    PData::PDataFact,
    dₙ::Vector,
    gt::Vector;
    λfact2::Bool = false,
)
    # printstyled("on est dans Shamanskii \n", color = :bold)
    nlp_at_x = nlp_stop.current_state
    x = nlp_at_x.x
    n = length(x)
    xdemi = x + dₙ
    gtemp = grad(nlp_stop.pb, xdemi)

    ϵ2 = sqrt(eps(eltype(gt)))
    Γ = max.(abs.(PData.Δ), ϵ2)

    if (PData.λ == 0) && !(λfact2)
        d̃ = PData.L \ gtemp[PData.pp]
        d̂ = PData.L' \ (PData.Q * (PData.Q' * d̃ ./ Γ))
        d = -d̂[invperm(PData.pp)]
    else
        # println("on est ici 🐵↓")
        # @show (nlp_at_x.Hx + PData.λ * Matrix(1.0I, n, n))
        L, D, pp, ρ, ncomp = ldlt_symm((nlp_at_x.Hx + PData.λ * Matrix(1.0I, n, n)))
        # println("On a L")
        D = Hermitian(D)
        X = eigen(D)
        Δ = X.values
        Q = X.vectors
        d̃ = L \ gtemp[pp]
        d̂ = L' \ (Q * (Q' * d̃ ./ Γ))
        d = -d̂[invperm(pp)]
    end

    dₕₒ = dₙ .+ d

    return dₕₒ
end
