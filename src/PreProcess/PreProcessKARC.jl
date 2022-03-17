function preprocessKARC(PData::PDataKARC, Hop, g, params::Tparams, calls, max_calls) #where T
    ζ = params.ζ
    nshifts = params.nshifts
    shifts = params.shifts

    n = length(g)
    gNorm2 = norm(g) # BLAS.nrm2(n, g, 1)
    precision = max(1e-12, min(0.5, (gNorm2^ζ)))
    ϵ = 1e-12#sqrt(eps()) # * 100.0
    cgtol = max(ϵ, min(0.09, 0.01 * norm(g)^(1.0 + ζ)))

    nshifts = length(shifts)
    solver = PData.solver
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

    PData.d .= g
    PData.λ = -1.0
    PData.ζ = ζ
    PData.indmin = 0
    PData.positives .= solver.converged
    for i=1:nshifts
        PData.xShift[i] .= solver.x[i]
        PData.norm_dirs[i] = norm(solver.x[i])
    end
    PData.shifts .= shifts
    PData.nshifts = nshifts
    PData.OK = true

    return PData # PDataKARC(g, -1.0, ζ, 0, positives, xShift, shifts, nshifts, Ndirs, true)
end
