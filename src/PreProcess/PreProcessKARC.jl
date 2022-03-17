function preprocessKARC(Hop, g, params::Tparams, calls, max_calls) #where T
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

    return PDataK(g, -1.0, ζ, 0, positives, xShift, shifts, nshifts, Ndirs, true)
end

function decreaseKARC(X::PDataK, α::Float64, TR::TrustRegion)
    X.indmin += 1
    p_imin = X.positives[X.indmin]
    α2 = max(X.norm_dirs[p_imin] / X.shifts[p_imin], eps())

    targetα = α * TR.decrease_factor
    @show X.positives
    @show X.indmin
    @show α, targetα, α2
    @show X.norm_dirs
    @show X.shifts
    @show X.norm_dirs ./ X.shifts
    @show X.xShift

    # fix α to its "ideal" value to satisfy αλ=||d||
    # while ensuring α decreases enough
    while α2 > targetα && p_imin < length(X.positives)
        X.indmin += 1
        p_imin = X.positives[X.indmin]
        α2 = max(X.norm_dirs[p_imin] / X.shifts[p_imin], eps())
    end

    if p_imin == length(X.positives)
        @warn "PreProcessKTR failure: α2=$α2"
    end

    X.d = X.xShift[p_imin]
    X.λ = X.shifts[p_imin]

    return α2
end
