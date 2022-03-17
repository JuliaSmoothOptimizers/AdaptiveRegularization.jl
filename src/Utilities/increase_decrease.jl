# default increase and decrease functions.
function decreaseBase(α::T, TR::TrustRegion) where {T}
    return α * TR.decrease_factor
end

function decrease(X::TPData, α::T, TR::TrustRegion) where {T}
    return decreaseBase(α, TR)
end

function decrease(X::PDataFact, α::T, TR::TrustRegion) where {T}
    X.success = false
    return decreaseBase(α, TR)
end

function decrease(X::PDataKARC, α::Float64, TR::TrustRegion)
    X.indmin += 1
    p_imin = X.positives[X.indmin]
    α2 = max(X.norm_dirs[p_imin] / X.shifts[p_imin], eps())

    targetα = α * TR.decrease_factor

    # fix α to its "ideal" value to satisfy αλ=||d||
    # while ensuring α decreases enough
    while α2 > targetα && p_imin < length(X.positives)
        X.indmin += 1
        p_imin = X.positives[X.indmin]
        α2 = max(X.norm_dirs[p_imin] / X.shifts[p_imin], eps())
    end

    if p_imin == length(X.positives)
        @warn "PreProcessKARC failure: α2=$α2"
    end

    X.d = X.xShift[p_imin]
    X.λ = X.shifts[p_imin]

    return α2
end

function decrease(X::PDataKTR, α::Float64, TR::TrustRegion)
    X.indmin += 1
    p_imin = X.positives[X.indmin]
    α2 = X.norm_dirs[p_imin]

    # fix α to its "ideal" value to satisfy α=||d||
    # while ensuring α decreases enough
    targetα = α * TR.decrease_factor

    while α2 > targetα && p_imin < length(X.positives)
        X.indmin += 1
        p_imin = X.positives[X.indmin]
        α2 = X.norm_dirs[p_imin]
    end

    if p_imin == length(X.positives)
        @warn "PreProcessKTR failure: α2=$α2"
    end

    X.d = X.xShift[p_imin]
    X.λ = X.shifts[p_imin]

    return α2
end

function increase(α::T, TR::TrustRegion) where {T}
    return min(α * TR.increase_factor, TR.max_α)
end

function increase(X::TPData, α::T, TR::TrustRegion) where {T}
    return increase(α, TR)
end
