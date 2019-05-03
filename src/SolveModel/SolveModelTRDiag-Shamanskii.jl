function solve_modelTRDiag_Shamanskii(PData :: PDataFact, δ:: Float64)
    # Solve the TR subproblem once diagonalized into Δ using the norm |Δ|
    # Setup the problem
    M = fill(1.0, size(PData.Δ))
    ϵ = 1.0e-10 * (1.0 + PData.λ)

    if PData.success # take care of eventual hard case and Newton's direction interior (λ = 0)
        # (PData.Δ + PData.λ * M) ⪰ 0
        λ = max(ϵ, PData.λ + ϵ) # to make sure (PData.Δ+λ*M) ≻ 0


        d̃₁ = -(PData.Δ .+ λ * M) .\ PData.g̃₁ # Ajouter Shamanskii ici!
        g̃₂ =  PData.Q' * grad(PData.f, Pdata.x + d̃₁)
        d̃₂ = -(PData.Δ .+ λ * M) .\ PData.g̃₂
        normd̃₁ = sqrt(d̃₂⋅d̃₂)
        if normd̃₁ < δ
            if PData.λ == 0.0 # Newton's direction
                λ = PData.λ
                ##println(" Newton's direction inside the region")
                #  d̃ is the Newton's direction, nothing more to do
            else              # hard case
                ##println(" hard case")
                bidon, i = findmin(PData.Δ)
                d̃₂[i] = 0.0; d̃₂[i] = - sign(PData.g̃­₂[i]) * sqrt(δ^2 - d̃₂⋅d̃₂)
            end
        else
            d̃₂,λ = solve_diagTR(λ, PData.Δ, PData.g̃, δ, ϵ)
        end
    else # hard case impossible, λ already > λ_min
        # Solve the subproblem (Δ + λ I) d̃ = -g̃₁ such that ||d̃|| = δ
        d̃₂, λ = solve_diagTR(PData.λ, PData.Δ, PData.g̃, δ, ϵ)
    end

    PData.λ = λ

    # Transform back d̃₁ into d
    d = AInv(PData,d̃₂)

    #try assert((PData.g̃₁ + 0.5*PData.Δ .* d̃₁)⋅d̃₁ <= 0.0)  catch  @bp  end

    return d, λ
end
