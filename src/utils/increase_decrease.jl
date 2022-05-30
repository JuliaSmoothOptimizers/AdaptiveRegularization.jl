function decrease(X::TPData, α::T, TR::TrustRegion) where {T}
    return α * TR.decrease_factor
end

function increase(::TPData, α::T, TR::TrustRegion) where {T}
    return min(α * TR.increase_factor, TR.max_α)
end

function decrease(X::PDataFact, α::T, TR::TrustRegion) where {T}
    X.success = false
    return α * TR.decrease_factor
end

# X.indmin is between 1 and length(positives)
# p_imin is between 1 and nshifts
function decrease(X::PDataKARC, α::Float64, TR::TrustRegion)
    positives = findall(X.positives)
    X.indmin += 1 # the step wasn't successful so we need to change something
    
    p_imin = positives[X.indmin]
    α2 = max(X.norm_dirs[p_imin] / X.shifts[p_imin], eps())

    targetα = α * TR.decrease_factor

    # fix α to its "ideal" value to satisfy αλ=||d||
    # while ensuring α decreases enough
    while α2 > targetα && X.indmin < length(positives)
        X.indmin += 1
        p_imin = positives[X.indmin]
        α2 = max(X.norm_dirs[p_imin] / X.shifts[p_imin], eps())
    end

    if X.indmin == length(positives) & (α2 > targetα) # p_imin == length(positives)
        @warn "PreProcessKARC failure: α2=$α2"
    end

    X.d = X.xShift[p_imin]
    X.λ = X.shifts[p_imin]

    return α2
end

function decrease(X::PDataKTR, α::Float64, TR::TrustRegion)
    X.indmin += 1
    positives = findall(X.positives)
    p_imin = positives[X.indmin]
    α2 = X.norm_dirs[p_imin]

    # fix α to its "ideal" value to satisfy α=||d||
    # while ensuring α decreases enough
    targetα = α * TR.decrease_factor

    while α2 > targetα && p_imin < length(positives)
        X.indmin += 1
        p_imin = positives[X.indmin]
        α2 = X.norm_dirs[p_imin]
    end

    if p_imin == length(positives)
        @warn "PreProcessKTR failure: α2=$α2"
    end

    X.d = X.xShift[p_imin]
    X.λ = X.shifts[p_imin]

    return α2
end
