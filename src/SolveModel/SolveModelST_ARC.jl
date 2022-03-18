function solve_modelST_ARC(H, g, nlp_stop, PData::PDataST, δ::T) where {T}
    # Steihaug-Toint adapted to ARC
    ϵ = sqrt(eps(T)) * 100.0
    n = length(g)
    cgtol = max(ϵ, min(0.7, 0.01 * norm(g)^(1.0 + PData.ζ)))

    (d, cg_stats) = cgARC(
        H,
        -g,
        atol = cgtol,
        rtol = 0.0,
        regulα = δ,
        itmax = max(2 * n, 50),
        verbose = false,
    )

    PData.d .= d
    return PData.d, PData.λ
end
