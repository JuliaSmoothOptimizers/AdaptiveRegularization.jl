function solve_modelARCDiagAbs(PData :: PDataFact, α:: Float64)
    # Solve the ARC subproblem once diagonalized into Δ using the norm |Δ|
    
    # Setup the subproblem to recover identity regularizing matrix
    # Should this be in pre-process so that PData.λ is indeed used as updated?
    # No No! Use PData.success to initialize PData.λ when true, and use it otherwise
    #
    # TO TEST to confirm it avoids ascent directions
    #
    M=ones(PData.Δ)
    n_g = norm(PData.g̃)
    ϵ =  1.0e-10 
    Γ = max(abs(PData.Δ),ϵ)

    S = sqrt(Γ)
    
    g_s = S .\ PData.g̃
    Δ = S .\ PData.Δ ./ S
    λmin = 0.0 #  Test the following
    if PData.success # ensure 
        l_m, = findmin(Δ)
        λ = max(-l_m,0.0)
        λmin = max(ϵ,  λ + ϵ * (1.0 + λ)) 
    else
        λmin = PData.λ
    end
    d_s = -(Δ+λmin*M) .\ g_s
    seuil_s = norm(d_s)/α


    # Solve the subproblem (Δ + λ I) d_s = -g_s such that λ = ||d_s||/α
    d_s,λ = solve_diag(λmin,Δ,g_s,seuil_s,α,ϵ)

    # Transform back d_s into d
    d̃ = S .\ d_s
    d = TtildeInv(PData,d̃)  
    
    #try assert((PData.g̃ + 0.5*PData.Δ .* d̃)⋅d̃ <= 0.0)  catch  @bp  end
    #println("*******SolveModelDiagAbs:  PData.g̃⋅d̃ = $(PData.g̃⋅d̃), 0.5 d̃'PData.Δd̃ = $(0.5*(PData.Δ .* d̃)⋅d̃)")

    PData.λ = λ

    return d, λ
end
