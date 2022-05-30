function preprocess(PData::PDataSpectral, H, g, n1, n2, α)

    X = eigen(H)
    Δ = X.values
    V = X.vectors

    l_m = minimum(Δ)
    g̃ = V' * g
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
