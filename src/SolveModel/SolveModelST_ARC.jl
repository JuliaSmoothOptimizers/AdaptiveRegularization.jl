function solve_modelST_ARC(X :: PDataST, δ:: Float64)
# Steihaug-Toint adapted to ARC
    e=1e-6
    n = length(X.g)
    cgtol = max(e, min(0.7, 0.01 * norm(X.g)^(1.0 + X.τ)))

    (d, cg_stats) = cgARC(X.H, -X.g,
                     atol=cgtol, rtol=0.0,
                     regulα=δ,
                     itmax=max(2 * n, 50),
                     verbose=false)
    λ = 0.0  #  dummy for this variant
    return d,λ
end
