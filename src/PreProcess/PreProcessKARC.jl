function preprocessKARC(Hop, g, params::Tparams, calls, max_calls) #where T
    τ = params.τ
    ζ = params.ζ
    nshifts = params.nshifts
    shifts = params.shifts

    n = length(g)
    gNorm2 = BLAS.nrm2(n, g, 1)
    precision = max(1e-12, min(0.5, (gNorm2^ζ)))
    ϵ = 1e-12#sqrt(eps()) # * 100.0
    cgtol = max(ϵ, min(0.09, 0.01 * norm(g)^(1.0 + ζ)))

    nshifts = length(shifts)
    solver = CgLanczosShiftSolver(Hop, -g, nshifts)
    cg_lanczos!(
        solver,
        Hop,
        -g,
        shifts,
        itmax = min(max_calls - sum(calls), 2 * n),
        atol = 1.0e-8, # cgtol
        rtol = precision, # ϵ
        verbose = 0,
        check_curvature = true,
    )
    xShift = solver.x
    positives = findall(.!solver.not_cv)
    Ndirs = [norm(dx) for dx in xShift]

    return PDataK(g, -1.0, ζ, τ, 0, positives, xShift, shifts, nshifts, Ndirs, true)
end

function decreaseKARC(X::PDataK, α::Float64, TR::TrustRegion)
    X.indmin += 1
    p_imin = X.positives[X.indmin]
    α2 = X.norm_dirs[p_imin]^X.τ / X.shifts[p_imin]

    targetα = α * TR.decrease_factor

    # fix α to its "ideal" value to satisfy αλ=||d||^τ
    # while ensuring α decreases enough
    while α2 > targetα && p_imin < length(X.positives)
        X.indmin += 1
        p_imin = X.positives[X.indmin]
        α2 = X.norm_dirs[p_imin]^X.τ / X.shifts[p_imin]
    end

    if p_imin == length(X.positives)
        @warn "PreProcessKTR failure no α2 found"
    end

    X.d = X.xShift[p_imin]
    X.λ = X.shifts[p_imin]

    return α2
end
