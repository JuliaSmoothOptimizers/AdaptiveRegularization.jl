function solve_modelST_TR_Sham_BFGS(nlp_stop, X :: PDataST, δ:: T) where T
# cas particulier Steihaug-Toint
    ϵ = sqrt(eps(T)) # * 100.0
    n = length(X.g)
    cgtol = max(ϵ, min(0.7, 0.01 * norm(X.g)^(1.0 + X.τ)))

    (dₙ, cg_stats_newt) = cg(X.H, -X.g, atol = cgtol, rtol = 1e-08,
                             radius = δ,  itmax = max(2 * n, 50),
                             verbose = false)



    λ = 0.0  #  dummy for this variant

    if cg_stats_newt.status == "solution good enough given atol and rtol"
        nlp_at_x = nlp_stop.current_state
        αₖ = 1.0
        sₖ = αₖ * dₙ
        xdemi = nlp_at_x.x + sₖ
        ∇fₖ₊₁ = grad(nlp_stop.pb, xdemi)
        yₖ = ∇fₖ₊₁ - nlp_at_x.gx

        Bₖ = X.H
        sBs = 1.0 / (sₖ' * Bₖ * sₖ)
        Yₖ = yₖ * yₖ'
        yts = yₖ' * sₖ
        Yₖ = Yₖ ./ yts
        Bₖ₊₁ = Bₖ + Yₖ + (Bₖ * sₖ * sₖ' * Bₖ') * sBs

        (dₛ, cg_stats) = cg(Bₖ₊₁, -∇fₖ₊₁, atol = cgtol, rtol = 1e-08,
                            radius = δ,  itmax = max(2 * n, 50),
                            verbose = false)

        Hxdₛ = X.H * dₛ
        if (norm(dₛ) < 2.0 * δ) && ((-(nlp_at_x.gx + 0.5 * Hxdₛ)⋅dₛ) > 0.0)
            return dₛ, xdemi, λ
        end
    end

    return dₙ, NaN * rand(length((dₙ))), λ
end
