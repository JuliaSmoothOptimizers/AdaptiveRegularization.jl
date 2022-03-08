export preprocessKTR, decreaseKTR
function preprocessKTR(Hop, g, params::Tparams, calls, max_calls)
    ζ = params.ζ
    #τ = params.τ
    nshifts = params.nshifts
    shifts = params.shifts

    n = length(g)
    gNorm2 = BLAS.nrm2(n, g, 1)
    precision =  max(1e-12,min(0.5,(gNorm2^ζ)))
    (xShift, stats) = cg_lanczos(Hop,
                                           -g,
                                           shifts,
                                           itmax=min(max_calls-sum(calls),2*n),
                                           #τ = τ,
                                           atol = 1.0e-8,
                                           rtol = precision,
                                           verbose=0,
                                           check_curvature=true,
                                           history = true)

    # positives = collect(findfirst(!, stats.flagged):length(stats.flagged))
    ϵ = 1.0e-8 + norm(g) * precision
    tmp = map(x -> x[end] <= ϵ, stats.residuals)
    positives = findall(tmp)
    if stats.solved
      @warn "We need to fix `stats.flagged`"
      # @show tmp, positives
    end

    success = false
    good_grad = false
    if VERSION >= v"1.1.0"
        xShift = hcat(xShift...)
    end
    dirs = [ (xShift[:,i]) for i = 1 : nshifts ];
    Ndirs = map(norm, dirs);

    d = g # bidon

    return  PDataK(d,-1.0, ζ, 0.0, 0, positives,xShift,shifts,nshifts,Ndirs,true)
end

function decreaseKTR(X :: PDataK, α:: Float64, TR:: TrustRegion)
    X.indmin += 1
    p_imin = X.positives[X.indmin]
    α2 = X.norm_dirs[p_imin]

    # fix α to its "ideal" value to satisfy α=||d||
    # while ensuring α decreases enough
    targetα = α*TR.decrease_factor

    while α2 > targetα
        X.indmin += 1
        p_imin = X.positives[X.indmin]
        α2 = X.norm_dirs[p_imin]
    end

    X.d = X.xShift[:,p_imin]
    X.λ = X.shifts[p_imin]

    return α2
end
