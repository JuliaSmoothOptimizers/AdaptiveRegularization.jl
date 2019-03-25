function solve_modelARCDiagAbs(PData :: PDataFact, α:: Float64)
    printstyled("on est dans solve_modelARCDiagAbs \n", color = :blue)
    # Solve the ARC subproblem once diagonalized into Δ using the norm |Δ|

    # Setup the subproblem to recover identity regularizing matrix
    # Should this be in pre-process so that PData.λ is indeed used as updated?
    # No No! Use PData.success to initialize PData.λ when true, and use it otherwise
    #
    # TO TEST to confirm it avoids ascent directions
    #
    ḡ = PData.g̃
    printstyled("on a ḡ = $ḡ\n", color = :yellow)
    n_g = norm(ḡ)
    printstyled("on a n_g = $n_g\n", color = :yellow)
    ϵ =  1.0e-10
    printstyled("on a ϵ = $ϵ \n", color = :yellow)
    Γ2 = max(abs(PData.Δ),ϵ)
    printstyled("on a Γ2 = $Γ2 \n", color = :yellow)
    Γ = sqrt(Γ2)
    printstyled("on a Γ = $Γ \n", color = :yellow)
    Δ̄ = Γ .\ PData.Δ ./ Γ
    printstyled("on a Δ̄  = $Δ̄  \n", color = :yellow)
    Δ = PData.Δ
    printstyled("on a Δ = $Δ \n", color = :yellow)
    λmin = 0.0 #  Test the following
    printstyled("λmin = $λmin \n", color = :yellow)
    if PData.success # ensure
        l_m, = findmin(Δ̄)
        λ = max(-l_m,0.0)
        λmin = max(ϵ,  λ + ϵ * (1.0 + λ))
    else
        λmin = PData.λ
    end
    d̄ = -(Δ + λmin*Γ2) .\ ḡ
    seuil_bar = norm(d̄)/α

    printstyled("on est a seuil_bar = $seuil_bar \n", color = :yellow)
    # Solve the subproblem (Δ̄ + λ I) d̄ = -ḡ such that λ = ||d̄||/α
    d̄,λ = solve_diag(λmin,Δ,ḡ,seuil_bar,α,ϵ, M=Γ2)
    # Transform back d_s into d
    d̃ = d̄
    #d̃ = Γ .\ d̄
    d = AInv(PData,d̃)

    #try assert((PData.g̃ + 0.5*PData.Δ .* d̃)⋅d̃ <= 0.0)  catch  @bp  end
    #println("*******SolveModelDiagAbs:  PData.g̃⋅d̃ = $(PData.g̃⋅d̃), 0.5 d̃'PData.Δd̃ = $(0.5*(PData.Δ .* d̃)⋅d̃)")

    PData.λ = λ

    return d, λ
end
