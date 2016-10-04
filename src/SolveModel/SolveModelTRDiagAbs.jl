function solve_modelTRDiagAbs(PData :: PDataFact, δ:: Float64)
    # Solve the TR subproblem once diagonalized into Δ using the norm |Δ|
    

    # Setup the subproblem to recover identity regularizing matrix
    # Should this be in pre-process so that PData.λ is indeed used as updated?
    #
    M=ones(PData.Δ)
    n_g = norm(PData.g̃)
    ϵ =  1.0e-10 # *n_g
    Γ = max(abs(PData.Δ),ϵ)

    S = sqrt(Γ)
    
    g_s = S .\ PData.g̃
    Δ = S .\ PData.Δ ./ S
    if PData.success    # initialize λmin, take care of hard case and Newton's direction interior
        l_m,i = findmin(Δ)
        λ = max(-l_m,0.0)
        λmin = max(ϵ,  λ + ϵ * (1.0 + λ)) #λ + ϵ #max(λ*(1+ϵ), ϵ)

        d_s = -(Δ+λmin*M) .\ g_s
        normd_s = sqrt(d_s⋅d_s)
        #if isnan(normd_s)  @bp end
        if normd_s < δ
            if λ != 0.0  # λ == 0 is Newton's direction, nothing to do...
                #println(" hard case")
                d_s[i] =0.0; d_s[i] = - sign(PData.g̃[i]) * sqrt(δ^2 - d_s⋅d_s)
            end
        else
            #println(" easy case normd_s = $normd_s")
            d_s,λ = solve_diagTR(λmin,Δ,g_s,δ,ϵ)
        end
    else # hard case impossible, λ already > λ_min
        # Solve the subproblem (Δ + λ I) d_s = -g_s such that δ = ||d_s||
        #println(" retry for success case")
        d_s,λ = solve_diagTR(PData.λ,Δ,g_s,δ,ϵ)
    end

    PData.λ = λ

    # Transform back d_s into d
    d̃ = S .\ d_s
    d = TtildeInv(PData,d̃)  
    #try assert((PData.g̃ + 0.5*PData.Δ .* d̃)⋅d̃ <= 0.0)  catch  @bp  end
#println("*******SolveModelDiagAbs:  PData.g̃⋅d̃ = $(PData.g̃⋅d̃), 0.5 d̃'PData.Δd̃ = $(0.5*(PData.Δ .* d̃)⋅d̃)")
    return d, λ
end