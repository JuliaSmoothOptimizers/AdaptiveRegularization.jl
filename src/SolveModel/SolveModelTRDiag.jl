export solve_modelTRDiag

function solve_modelTRDiag(nlp_stop, PData :: PDataFact, δ:: T) where T
    # Solve the TR subproblem once diagonalized into Δ using the norm |Δ|
    # Setup the problem
    # println("On est ici!")
    M = T.(fill(1.0, size(PData.Δ)))
    ϵ = sqrt(eps(T)) / 100.0
    ϵ2 = T.(ϵ * (1.0 + PData.λ))


    if PData.success # take care of eventual hard case and Newton's direction interior (λ = 0)
        # (PData.Δ + PData.λ * M) ⪰ 0
        λ = max(ϵ2, PData.λ + ϵ2) # to make sure (PData.Δ + λ * M) ≻ 0

        d̃ = -(PData.Δ .+ λ * M) .\ PData.g̃ # Ajouter Shamanskii ici!
        normd̃ = sqrt(d̃⋅d̃)
        if normd̃ < δ
            if PData.λ == 0.0 # Newton's direction
                λ = PData.λ
                # println(" Newton's direction inside the region")
                #  d̃ is the Newton's direction, nothing more to do
            else              # hard case
                ##println(" hard case")
                bidon, i = findmin(PData.Δ)
                d̃[i] = 0.0; d̃[i] = - sign(PData.g̃[i]) * sqrt(δ^2 - d̃⋅d̃)
            end
        else
            d̃,λ = solve_diagTR(λ, PData.Δ, PData.g̃, δ, ϵ)
        end
    else # hard case impossible, λ already > λ_min
        # println("on est dans le cas !(PData.success)")
        # Solve the subproblem (Δ + λ I) d̃ = -g̃ such that ||d̃|| = δ
        d̃, λ = solve_diagTR(PData.λ, PData.Δ, PData.g̃, δ, ϵ)
    end

    PData.λ = λ

    # Transform back d̃ into d
    d = AInv(PData, d̃)

    #try assert((PData.g̃ + 0.5*PData.Δ .* d̃)⋅d̃ <= 0.0)  catch  @bp  end
    return d, NaN * ones(size(d)), λ
end
