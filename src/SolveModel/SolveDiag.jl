function solve_diag(λ,Δ,g̃,seuil,α,ϵ)

    M = ones(Δ)
    λin = λ
    λtry = seuil
    dtry = -(Δ+λtry*M) .\ g̃
    normdtry = sqrt(dtry⋅dtry)
    seuiltry = normdtry/α
    d̃ = dtry
    #biss_iter = 0
    #println(" Biss $biss_iter : λ = $λtry  seuil = $seuiltry ")

    # secant like heuristic to help Newton process by providing a good starting point
    iter_bis= 0
    if λin < seuil
        while λtry > seuiltry
            λtry = λ + (λtry - λ)/max((seuil-seuiltry),2.0)
            dtry = -(Δ+λtry*M) .\ g̃
            normdtry = sqrt(dtry⋅dtry)
            seuiltry = normdtry/α
        end
        λ = λtry
        d̃ = dtry
        normd̃ = normdtry
        seuil = seuiltry
        iter_bis += 1
    end

    
    # Newton iterations
    iter_nwt = 0
    while (abs(seuil - λ) > ϵ) && (iter_nwt<40)
        dotd̃ = (Δ+λ*M) .\ d̃
        Δλ = λ*(seuil - λ) / (seuil + λ*(λ*(d̃ ⋅ dotd̃)/(d̃ ⋅ d̃)))
        λ =  max(λ + Δλ, λin+1.0e-10)
        d̃ = -(Δ+λ*M) .\ g̃
        normd̃ = sqrt(d̃⋅d̃)
        seuil = normd̃/α

        assert(λ>=0.0)
        iter_nwt += 1
    end

    #try assert((g̃ + 0.5*Δ.*d̃)⋅d̃ <= 0.0)  catch  @bp  end
    #println("*******SolveDiag:  g̃⋅d̃ = $(g̃⋅d̃), 0.5 d̃'Δd̃ = $(0.5*(Δ .* d̃)⋅d̃)")

    return d̃,λ
end
