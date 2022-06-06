function solve_modelKARC(H, g, nlp_stop, X::PDataKARC, α::T) where {T}
    # target value should be close to satisfy αλ=||d||
    start = findfirst(X.positives)
    if isnothing(start)
        start = 1
    end
    positives = if VERSION < v"1.7.0"
        collect(start:length(X.positives))
    else
        start:length(X.positives)
    end
    target = ((abs(α * X.shifts[i] - X.norm_dirs[i])) for i in positives)

    # pick the closest shift to the target within positive definite H+λI
    indmin = argmin(target)
    X.indmin = start + indmin - 1
    p_imin = X.indmin
    X.d .= X.xShift[p_imin]
    X.λ = X.shifts[p_imin]

    return X.d, X.λ
end
