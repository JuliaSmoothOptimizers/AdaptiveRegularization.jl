function solve_modelST_TR(H, g, nlp_stop, PData::PDataST, δ::T; cgtol::T = 0.1) where {T}
    # cas particulier Steihaug-Toint
    # ϵ = sqrt(eps(T)) # * 100.0 # old
    # cgtol = max(ϵ, min(cgtol, 9 * cgtol / 10, 0.01 * norm(g)^(1.0 + PData.ζ))) # old

    ζ, ξ, maxtol, mintol = PData.ζ, PData.ξ, PData.maxtol, PData.mintol
    n = length(g)
    gNorm2 = norm(g)
    # precision = max(1e-12, min(0.5, (gNorm2^ζ)))
    # Tolerance used in Assumption 2.6b in the paper ( ξ > 0, 0 < ζ ≤ 1 )
    cgatol = min(maxtol, ξ * gNorm2, ξ * gNorm2^(1 + ζ))
    cgatol = max(mintol, cgatol) # add some feasible limit
    cgrtol = eps(eltype(g))

    solver = PData.solver
    cg!(
        solver,
        H,
        -g,
        atol = cgatol,
        rtol = cgrtol,
        radius = δ,
        itmax = max(2 * n, 50),
        verbose = 0,
    )

    PData.d .= solver.x
    PData.OK = solver.stats.solved

    return PData.d, PData.λ
end
