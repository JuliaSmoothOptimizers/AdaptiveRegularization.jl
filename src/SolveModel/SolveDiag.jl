function solve_diag(λ,Δ,g̃,seuil,α,ϵ; M = [0.0])

    if M==[0.0]  M = ones(Δ) ; end
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
        while λtry > seuiltry && (iter_bis<10)
            λtry = λ + (λtry - λ)/max((seuil-seuiltry),2.0)
            dtry = -(Δ+λtry*M) .\ g̃
            normdtry = sqrt(dtry⋅dtry)
            seuiltry = normdtry/α
            iter_bis += 1
#println(" λ computation, iter_bis : $iter_bis  λtry  $λtry")
        end
        λ = λtry
        d̃ = dtry
        normd̃ = normdtry
        seuil = seuiltry
    end

    
    # Newton iterations
    iter_nwt = 0
    while (abs(seuil - λ)/max(seuil,λ) > ϵ) && (iter_nwt<40)
        dotd̃ = (Δ+λ*M) .\ (M .* d̃)
        Δλ = λ*(seuil - λ) / (seuil + λ*(λ*(d̃ ⋅ dotd̃)/(d̃ ⋅ d̃)))
        λ =  max(λ + Δλ, λin+1.0e-10)
        d̃ = -(Δ+λ*M) .\ g̃
        normd̃ = sqrt(d̃⋅d̃)
        seuil = normd̃/α
#println(" λ computation,  iter_nwt : $iter_nwt  λ  $λ")

        assert(λ>=0.0)
        iter_nwt += 1
    end
#println(" λ computation, iter_bis : $iter_bis  iter_nwt : $iter_nwt  λ  $λ")
    #try assert((g̃ + 0.5*Δ.*d̃)⋅d̃ <= 0.0)  catch  @bp  end
    #println("*******SolveDiag:  g̃⋅d̃ = $(g̃⋅d̃), 0.5 d̃'Δd̃ = $(0.5*(Δ .* d̃)⋅d̃)")

    return d̃,λ
end
