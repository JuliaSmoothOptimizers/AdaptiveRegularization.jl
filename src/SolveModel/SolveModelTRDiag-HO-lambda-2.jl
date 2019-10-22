export solve_modelTRDiag_HO_λ_2

"""
If the Newton direction is accepted and the high order correction lies within
a bigger trust region then we use the high order correction.
"""
function solve_modelTRDiag_HO_λ_2(nlp_stop, PData :: PDataFact, δ:: T; ho_correction :: Symbol = :Shamanskii, fact = 2.0, λfact = Inf) where T
    # Solve the TR subproblem once diagonalized into Δ using the norm |Δ|
    # Setup the problem
    # printstyled("On est dans solve_modelTRDiag_HO_λ ↓ \n", color = :red)
    nlp_at_x = nlp_stop.current_state
    M = fill(T.(1.0), size(PData.Δ))
    ϵ = sqrt(eps(T)) / T(100.0)
    ϵ2 = T.(ϵ * (T(1.0) + PData.λ))
    global dTR = nothing
    global dHO = nothing

    if PData.success # take care of eventual hard case and Newton's direction interior (λ = 0)
        # printstyled("on a PData.succes = $(PData.success) \n", color = :red)
        # (PData.Δ + PData.λ * M) ⪰ 0
        λ = max(ϵ2, PData.λ + ϵ2) # to make sure (PData.Δ + λ * M) ≻ 0

        d̃ = -(PData.Δ .+ λ * M) .\ PData.g̃
        normd̃ = sqrt(d̃⋅d̃)
        normg̃ = sqrt(PData.g̃⋅PData.g̃)
        if normd̃ < δ
            if PData.λ == 0.0 # Newton's direction
                # println(" Newton's direction inside the region")
                λ = PData.λ
                dN = AInv(PData, d̃)
                dHO = eval(ho_correction)(nlp_stop, PData, dN, PData.g̃)
                if (norm(dHO) < 2.0 .* δ) && ((-(nlp_at_x.gx + 0.5 * nlp_at_x.Hx * dHO)⋅dHO) > 0.0) # && (norm(grad(nlp_stop.pb, nlp_at_x.x + dHO)) < norm(nlp_at_x.gx))
                    # printstyled("on prend dHO 🐣\n", color = :green)
                    return dHO, dHO, λ
                else
                    return dN, dHO, λ
                end
                # return dN, dHO, xdemi, λ
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
        # printstyled("hard case impossible \n", color = :yellow)
        d̃, λ = solve_diagTR(PData.λ, PData.Δ, PData.g̃, δ, ϵ)
    end

    PData.λ = λ

    # Transform back d̃ into d
    d = AInv(PData, d̃)
    # printstyled("On a d après AInv \n", color = :red)
    rₙ = (-(nlp_at_x.gx + 0.5 * nlp_at_x.Hx * d)⋅d)
    if PData.λ < λfact
        # println("ici ! 🍆")
        # @show λfact
        # @show ho_correction
        dHO = eval(ho_correction)(nlp_stop, PData, d, PData.g̃)
        rₕₒ = (-(nlp_at_x.gx + 0.5 * nlp_at_x.Hx * dHO)⋅dHO)
        if (norm(dHO) < 2.0 .* δ) && (rₕₒ > rₙ)
            # printstyled("on prend dHO 🐣\n", color = :green)
            return dHO, dHO, PData.λ
        else
            # printstyled("on prend d 🐲 \n", color = :green)
            return d, dHO, PData.λ
        end
    end
    return d, NaN * rand(length(d)), PData.λ
    # return d, PData.λ
end
