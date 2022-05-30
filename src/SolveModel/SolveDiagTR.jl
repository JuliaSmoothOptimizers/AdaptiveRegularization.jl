function solve_diagTR(λ, Δ, g̃, δ, ϵ; M = [0.0])
    # λ underestimates the required value λstar such that ||d̃(λstar)|| = δ
    # where d̃(λ) =  -(Δ + λ * M) .\ g̃
    T = eltype(δ)
    M = T.(M)

    if M == [0.0]
        M = fill(T.(1.0), size(Δ)[1])
    end
    @assert (minimum(M) > 0.0)
    λin = T.(λ)
    d̃ = -(Δ .+ λ * M) .\ g̃
    d̃d̃ = d̃ ⋅ d̃
    normd̃ = sqrt(d̃d̃)

    tol1 = sqrt(eps(T))
    tolerance = tol1 * T(100.0)

    # Newton iterations
    iter_nwt = 0
    while (normd̃ >= (δ + tolerance * normd̃)) && (iter_nwt < 40)
        dotd̃ = (Δ .+ λ * M) .\ (M .* d̃)
        Δλ = ((normd̃ - δ) / δ) * (d̃d̃ / (d̃ ⋅ dotd̃))
        λeps = tol1 / T(100.0)
        λ = max(λ + Δλ, λin + λeps)
        d̃ = -(Δ .+ λ * M) .\ g̃
        d̃d̃ = d̃ ⋅ d̃
        normd̃ = sqrt(d̃d̃)
        @assert (λ >= 0.0)
        iter_nwt += 1
    end

    return d̃, λ
end
