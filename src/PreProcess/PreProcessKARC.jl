function preprocess(PData::PDataKARC, Hop, g, calls, max_calls, α)
    ζ, ξ, maxtol, mintol = PData.ζ, PData.ξ, PData.maxtol, PData.mintol
    nshifts = PData.nshifts
    shifts = PData.shifts

    n = length(g)
    gNorm2 = norm(g)
    # precision = max(1e-12, min(0.5, (gNorm2^ζ)))
    # Tolerance used in Assumption 2.6b in the paper ( ξ > 0, 0 < ζ ≤ 1 )
    cgatol = PData.cgatol(ζ, ξ, maxtol, mintol, gNorm2)
    cgrtol = PData.cgrtol(ζ, ξ, maxtol, mintol, gNorm2)

    nshifts = length(shifts)
    cb = (slv, A, b, shifts) -> begin
        ind = setdiff(1:length(shifts), findall(slv.not_cv))
        target = [norm(slv.x[i])/shifts[i] - α for i in ind] # the last one should be negative
        if length(ind) > 1
        if !isnothing(findfirst(target .> 0))
            slv.stats.user = true
        end
        end
    end
    solver = PData.solver
    cg_lanczos!(
        solver,
        Hop,
        -g,
        shifts,
        itmax = min(max_calls - sum(calls), 2 * n),
        atol = cgatol,
        rtol = cgrtol,
        verbose = 0,
        check_curvature = true,
        callback = cb,
    )

    PData.indmin = 0
    PData.positives .= solver.converged
    for i = 1:nshifts
        PData.xShift[i] .= solver.x[i]
        PData.norm_dirs[i] = norm(solver.x[i])
    end
    PData.shifts .= shifts
    PData.nshifts = nshifts
    PData.OK = sum(solver.converged) != 0 # at least one system was solved

    PData.solver.stats.user = false

    return PData
end
