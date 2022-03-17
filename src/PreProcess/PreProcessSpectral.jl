function preprocessSpectral(PData::PDataSpectral, H, g, n1, n2)
    #global V, D, rho, ncomp

    # Δ, V = eig(H)
    X = eigen(H)
    Δ = X.values
    V = X.vectors

    l_m, = findmin(Δ)
    g̃ = V' * g
    n_g = norm(g)
    ɛ = 1.0e-10
    λ = max(-l_m, 0.0)

    PData.V = V
    PData.Δ = Δ
    PData.g̃ = g̃
    PData.λ = λ
    PData.success = true
    PData.OK = true

    return PData
end


function AInv(X::PDataSpectral, d̃::Array{Float64,1})
    return X.V * d̃
end
