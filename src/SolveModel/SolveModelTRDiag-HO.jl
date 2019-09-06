export solve_modelTRDiag_HO

"""
If the Newton direction is accepted and the high order correction lies within
a bigger trust region then we use the high order correction.
"""
function solve_modelTRDiag_HO(nlp_stop, PData :: PDataFact, δ:: T; ho_correction :: Symbol = :Shamanskii, fact = 2.0) where T
    # Solve the TR subproblem once diagonalized into Δ using the norm |Δ|
    # Setup the problem
    # printstyled("On est dans solve_modelTRDiag_HO \n")
    nlp_at_x = nlp_stop.current_state
    M = fill(T.(1.0), size(PData.Δ))
    ϵ = sqrt(eps(T)) / T(100.0)
    ϵ2 = T.(ϵ * (T(1.0) + PData.λ))

    # @show PData.g̃e

    if PData.success # take care of eventual hard case and Newton's direction interior (λ = 0)
        # (PData.Δ + PData.λ * M) ⪰ 0
        λ = max(ϵ2, PData.λ + ϵ2) # to make sure (PData.Δ + λ * M) ≻ 0
        # @show λ
        # @show PData.Δ
        # @show M

        d̃ = -(PData.Δ .+ λ * M) .\ PData.g̃
        # @show d̃
        normd̃ = sqrt(d̃⋅d̃)
        # @show normd̃
        normg̃ = sqrt(PData.g̃⋅PData.g̃)
        # @show normg̃
        if normd̃ < δ
            if PData.λ == 0.0 # Newton's direction
                λ = PData.λ
                # @show λ
                dN = AInv(PData, d̃)
                # @show ho_correction
                # @show dN
                # @show PData.g̃
                dtemp, xdemi = eval(ho_correction)(nlp_stop, PData, dN, PData.g̃)
                dHO = dtemp
                # if (norm(dHO) < fact * δ) && ((-(nlp_at_x.gx + 0.5 * nlp_at_x.Hx * dHO)⋅dHO) > 0.0)
                if ((-(nlp_at_x.gx + 0.5 * nlp_at_x.Hx * dHO)⋅dHO) > 0.0)
                    return dHO, xdemi, λ
                else
                    return dN, xdemi, λ
                end
                return dN, dHO, xdemi, λ
                # println(" Newton's direction inside the region")
                #  d̃ is the Newton's direction, nothing more to do
            else              # hard case
                # println(" hard case")
                bidon, i = findmin(PData.Δ)
                d̃[i] = T(0.0); d̃[i] = - sign(PData.g̃[i]) * sqrt(δ^2 - d̃⋅d̃)
            end
        else
            # println("success pas hard case")
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
    return d, NaN * rand(length(d)), λ
end
