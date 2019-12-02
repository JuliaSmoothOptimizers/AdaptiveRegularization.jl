function preprocessKARC(Hop, g, params::Tparams, calls, max_calls) #where T
    τ = params.τ
    ζ = params.ζ
    nshifts = params.nshifts
    shifts = params.shifts

    n = length(g)
    gNorm2 = BLAS.nrm2(n, g, 1)
    precision =  max(1e-12,min(0.5,(gNorm2^ζ)))

    #ϵ = sqrt(eps(T)) # * 100.0
    ϵ = 1e-12#sqrt(eps()) # * 100.0
    cgtol = max(ϵ, min(0.09, 0.01 * norm(g)^(1.0 + ζ)))


    (xShift, stats) = cg_lanczos_shift_seq(Hop,
                                           -g,
                                           shifts,
                                           itmax=min(max_calls-sum(calls),2*n),
                                           #τ = τ,
                                           atol = 1.0e-8,
                                           rtol = precision,
                                           verbose=false,
                                           check_curvature=true)

#    (xShift, stats) = cg_lanczos_shift_seq(Hop,
#                                           -g,
#                                           shifts,
#                                           itmax=min(max_calls-sum(calls),2*n),
#                                           #τ = τ,
#                                           atol = cgtol,
#                                           rtol = ϵ,
#                                           verbose=false,
#                                           check_curvature=true)


    positives = collect(findfirst(!, stats.flagged):length(stats.flagged))

    success = false
    good_grad = false
    if VERSION >= v"1.1.0"
        xShift = hcat(xShift...)
    end
    dirs = [ (xShift[:,i]) for i = 1 : nshifts ];
    Ndirs = map(norm, dirs);

    d = g # bidon
    return  PDataK(d, -1.0, ζ, τ, 0 , positives, xShift,   shifts, nshifts, Ndirs, true)
end


function decreaseKARC(X :: PDataK, α:: Float64, TR:: TrustRegion)
    X.indmin += 1
    p_imin = X.positives[X.indmin]
    α2 = X.norm_dirs[p_imin]^X.τ/X.shifts[p_imin]

    targetα = α*TR.decrease_factor

    # fix α to its "ideal" value to satisfy αλ=||d||^τ
    # while ensuring α decreases enough
    while α2 > targetα
        X.indmin += 1
        p_imin = X.positives[X.indmin]
        α2 = X.norm_dirs[p_imin]^X.τ/X.shifts[p_imin]
    end

    X.d = X.xShift[:,p_imin]
    X.λ = X.shifts[p_imin]

    return α2
end
