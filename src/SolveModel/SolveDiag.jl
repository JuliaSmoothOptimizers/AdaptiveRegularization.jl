function solve_diag(λ, Δ, g̃, seuil, α::T, ϵ; M = [0.0]) where {T}
    if M == [0.0]
        M = ones(T, size(Δ))
    end
    λin = λ
    λtry = seuil
    dtry = -(Δ .+ λtry * M) .\ g̃
    normdtry = sqrt(dtry ⋅ dtry)
    seuiltry = normdtry / α
    d̃ = dtry

    # secant like heuristic to help Newton process by providing a good starting point
    global iter_bis = 0
    if λin < seuil
        while (λtry > seuiltry && (iter_bis < 10))
            λtry = λ + (λtry - λ) / max((seuil - seuiltry), T(2.0))
            dtry = -(Δ .+ λtry * M) .\ g̃
            normdtry = sqrt(dtry ⋅ dtry)
            seuiltry = normdtry / α
            global iter_bis += 1
        end
        λ = λtry
        d̃ = dtry
        normd̃ = normdtry
        seuil = seuiltry
    end


    # Newton iterations
    global iter_nwt = 0
    while (abs(seuil - λ) / max(seuil, λ) > ϵ) && (iter_nwt < 40)
        dotd̃ = (Δ .+ λ * M) .\ (M .* d̃)
        Δλ = λ * (seuil - λ) / (seuil + λ * (λ * (d̃ ⋅ dotd̃) / (d̃ ⋅ d̃)))
        λ = max(λ + Δλ, λin + T(1.0e-10))
        λ = max(λ + Δλ, λin + (sqrt(eps(T) ./ 100.0)))
        d̃ = -(Δ .+ λ * M) .\ g̃
        normd̃ = sqrt(d̃ ⋅ d̃)
        seuil = normd̃ / α


        @assert (λ >= 0.0)
        global iter_nwt += 1
    end
    return d̃, λ
end
