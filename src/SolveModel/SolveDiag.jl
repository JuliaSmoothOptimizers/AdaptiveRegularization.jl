function solve_diag(λ, Δ, g̃, seuil, α, ϵ; M = [0.0])
    # printstyled("on est dans solve_diag \n", color = :cyan)
    # printstyled("on a M = $M \n", color = :yellow)
    if M == [0.0]
        # M = ones(Δ) ;
        M = fill(1.0, size(Δ))
    end
    λin = λ
    # printstyled("on a λin = $λin \n", color = :yellow)
    λtry = seuil
    # printstyled("on a λtry = $λtry \n", color = :yellow)
    dtry = -(Δ .+ λtry * M) .\ g̃
    # printstyled("on a dtry = $dtry \n", color = :yellow)
    normdtry = sqrt(dtry⋅dtry)
    # printstyled("on a normdtry = $normdtry \n", color = :yellow)
    seuiltry = normdtry / α
    # printstyled("on a seuiltry = $seuiltry \n", color = :yellow)
    d̃ = dtry
    # printstyled("on a d̃ = $d̃ \n", color = :yellow)
    #biss_iter = 0
    #println(" Biss $biss_iter : λ = $λtry  seuil = $seuiltry ")

    # secant like heuristic to help Newton process by providing a good starting point
    global iter_bis = 0
    # printstyled("on a iter_bis = $iter_bis \n", color = :yellow)
    if λin < seuil
        # printstyled("on a λin < seuil = $(λin < seuil) \n", color = :yellow)
        while (λtry > seuiltry && (iter_bis < 10))
            # printstyled("on est dans le while \n", color = :bold)
            λtry = λ + (λtry - λ) / max((seuil - seuiltry), 2.0)
            dtry = -(Δ .+ λtry * M) .\ g̃
            normdtry = sqrt(dtry⋅dtry)
            seuiltry = normdtry / α
            global iter_bis += 1
            # printstyled("iter_bis = $iter_bis et λtry = $λtry et seuiltry = $seuiltry \n", color = :yellow)
            #println(" λ computation, iter_bis : $iter_bis  λtry  $λtry")
        end
        λ = λtry
        # printstyled("on a λ = $λ \n", color = :yellow)
        d̃ = dtry
        # printstyled("on a d̃ = $d̃ \n", color = :yellow)
        normd̃ = normdtry
        # printstyled("on a normd̃ = $normd̃ \n", color = :yellow)
        seuil = seuiltry
        # printstyled("on a seuil = $seuil \n", color = :yellow)
    end


    # Newton iterations
    global iter_nwt = 0
    # printstyled("on a iter_nwt = $iter_nwt \n", color = :yellow)
    while (abs(seuil - λ)/max(seuil,λ) > ϵ) && (iter_nwt < 40)
        # printstyled("On est dans le 2e while \n", color = :bold)
        # @show Δ
        # @show λ
        # @show M
        # # @show d̃
        dotd̃ = (Δ .+ λ * M) .\ (M .* d̃)
        # printstyled("on a dotd̃ = $dotd̃ \n", color = :yellow)
        Δλ = λ * (seuil - λ) / (seuil + λ * (λ * (d̃ ⋅ dotd̃) / (d̃ ⋅ d̃)))
        # printstyled("on a Δλ = $Δλ \n", color = :yellow)
        λ =  max(λ + Δλ, λin + 1.0e-10)
        # printstyled("on a λ = $λ \n", color = :yellow)
        d̃ = -(Δ .+ λ * M) .\ g̃
        # printstyled("on a d̃ = $d̃ \n", color = :yellow)
        normd̃ = sqrt(d̃⋅d̃)
        # printstyled("on a normd̃ = $normd̃ \n", color = :yellow)
        seuil = normd̃ / α
        # printstyled("on a α  = $α et seuil = $seuil \n", color = :yellow)

        # printstyled("avant @assert λ = $λ \n", color = :yellow)

        # @assert (λ >= 0.0)
        # printstyled("après @assert λ >= 0.0 = $(λ >= 0.0) \n", color = :yellow)
        global iter_nwt += 1
        # printstyled("on iter_nwt = $iter_nwt \n", color = :yellow)
    end
    #println(" λ computation, iter_bis : $iter_bis  iter_nwt : $iter_nwt  λ  $λ")
    #try assert((g̃ + 0.5*Δ.*d̃)⋅d̃ <= 0.0)  catch  @bp  end
    #println("*******SolveDiag:  g̃⋅d̃ = $(g̃⋅d̃), 0.5 d̃'Δd̃ = $(0.5*(Δ .* d̃)⋅d̃)")

    # printstyled("on a d̃ = $d̃ et λ = $λ dans solve_diag \n", color = :yellow)
    return d̃,λ
end
