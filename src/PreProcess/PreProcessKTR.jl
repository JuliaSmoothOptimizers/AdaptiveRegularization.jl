function preprocess(PData::PDataKTR, Hop, g, calls, max_calls)
    ζ = PData.ζ
    nshifts = PData.nshifts
    shifts = PData.shifts

    n = length(g)
    gNorm2 = norm(g)
    precision = max(1e-12, min(0.5, (gNorm2^ζ)))

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

    PData.indmin = 0
    PData.positives .= solver.converged
    for i=1:nshifts
        PData.xShift[i] .= solver.x[i]
        PData.norm_dirs[i] = norm(solver.x[i])
    end
    PData.shifts .= shifts
    PData.nshifts = nshifts
    PData.OK = sum(solver.converged) != 0 # at least one system was solved

    return PData
end
