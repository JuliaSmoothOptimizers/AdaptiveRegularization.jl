function solve_modelST_TR(nlp_stop, X :: PDataST, δ:: T) where T
# cas particulier Steihaug-Toint
    # e=1e-6
    ϵ = sqrt(eps(T)) * 100.0
    n = length(X.g)
    cgtol = max(e, min(0.7, 0.01 * norm(X.g)^(1.0 + X.τ)))

    (d, cg_stats) = cg(X.H, -X.g, atol = cgtol, rtol = 0.0,
                     radius = δ,  itmax = max(2 * n, 50),
                     verbose = false)

    λ = 0.0  #  dummy for this variant
    return d, NaN * rand(length((d))), λ
end
