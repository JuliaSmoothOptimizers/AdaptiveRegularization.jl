function solve_modelST_TR(H, g, nlp_stop, PData::PDataST, δ::T; cgtol::T = 0.1) where {T}
    # cas particulier Steihaug-Toint
    ϵ = sqrt(eps(T)) # * 100.0
    n = length(g)
    cgtol = max(ϵ, min(cgtol, 9 * cgtol / 10, 0.01 * norm(g)^(1.0 + PData.ζ)))

    solver = PData.solver
    cg!(
        solver,
        H,
        -g,
        atol = cgtol,
        rtol = ϵ,
        radius = δ,
        itmax = max(2 * n, 50),
        verbose = 0,
    )

    PData.d .= solver.x
    PData.OK = solver.stats.solved

    return PData.d, PData.λ
end
