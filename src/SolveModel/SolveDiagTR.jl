function solve_diagTR(λ, Δ, g̃, δ, ϵ; M = [0.0])
    # λ underestimates the required value λstar such that ||d̃(λstar)|| = δ
    # where d̃(λ) =  -(Δ + λ * M) .\ g̃
    # printstyled("dans solve_diagTR T(λ) = $(typeof(λ)) \n", color = :cyan)
    # printstyled("dans solve_diagTR T(Δ) = $(typeof(Δ)) \n", color = :cyan)
    # printstyled("dans solve_diagTR T(g̃) = $(typeof(g̃)) \n", color = :cyan)
    # printstyled("dans solve_diagTR T(δ) = $(typeof(δ)) \n", color = :cyan)
    # printstyled("dans solve_diagTR T(M) = $(eltype(M)) \n", color = :cyan)
    # @show M
    T = eltype(δ)
    # @show T
    M = T.(M)
    # @show M
    # printstyled("dans solve_diagTR T(M) = $(eltype(M)) \n", color = :cyan)

    if M == [0.0]
        # M = ones(Δ) ;
        M = fill(T.(1.0), size(Δ)[1]);
    end
    @assert (minimum(M) > 0.0)
    λin = T.(λ)
    # printstyled("dans solve_diagTR T(λin) = $(eltype(λin)) \n", color = :cyan)
    # @show Δ
    # @show λ
    # @show M
    # @show g̃
    d̃ = -(Δ .+ λ * M) .\ g̃
    # printstyled("dans solve_diagTR T(d̃) = $(eltype(d̃)) \n", color = :cyan)
    d̃d̃ = d̃⋅d̃
    # printstyled("dans solve_diagTR T(d̃d̃) = $(eltype(d̃d̃)) \n", color = :cyan)
    normd̃ = sqrt(d̃d̃)
    # printstyled("dans solve_diagTR T(normd̃) = $(eltype(normd̃)) \n", color = :cyan)

    # tolerance = 1e-06
    # @show sqrt(eps(T))
    tol1 = sqrt(eps(T))
    tolerance = tol1 * T(100.0)
    # printstyled("dans solve_diagTR T(tolerance) = $(eltype(tolerance)) \n", color = :cyan)
    # Newton iterations
    iter_nwt = 0
    #println(" Nwt $iter_nwt : λ = $λ  normd̃ = $normd̃  δ = $δ)")
    while (normd̃ >= (δ + tolerance * normd̃)) && (iter_nwt < 40)
        dotd̃ = (Δ .+ λ * M) .\ (M .* d̃)
        Δλ = ((normd̃-δ)/δ) * (d̃d̃/(d̃ ⋅ dotd̃))
        λeps = tol1 / T(100.0)
        λ =  max(λ + Δλ, λin + λeps)
        # printstyled("dans solve_diagTR dans le while λ = $λ \n", color = :cyan)
        # printstyled("dans solve_diagTR dans le while eltype(λ) = $(eltype(λ)) \n", color = :cyan)
        d̃ = -(Δ .+ λ * M) .\ g̃
        d̃d̃ = d̃⋅d̃
        normd̃ = sqrt(d̃d̃)
        @assert (λ >= 0.0)
        iter_nwt += 1
    end
    #println(" λ computation, iter_nwt : $iter_nwt  λ  $λ")
    # printstyled("dans solve_diagTR T(d̃) = $(eltype(d̃)) \n", color = :cyan)
    return d̃, λ
end
