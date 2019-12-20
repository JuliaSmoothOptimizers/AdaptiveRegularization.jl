function solve_modelST_ARC(nlp_stop, X :: PDataST, δ:: T) where T
# Steihaug-Toint adapted to ARC
    ϵ = sqrt(eps(T)) * 100.0
    n = length(X.g)
    cgtol = max(ϵ, min(0.7, 0.01 * norm(X.g)^(1.0 + X.ζ)))

    (d, cg_stats) = cgARC(X.H, -X.g,
                     atol = cgtol, rtol=0.0,
                     regulα = δ,
                     itmax = max(2 * n, 50),
                     verbose=false)

    λ = 0.0  #  dummy for this variant
    return d, NaN * rand(length(d)), λ
end
