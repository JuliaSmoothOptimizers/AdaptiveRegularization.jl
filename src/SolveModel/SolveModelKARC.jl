function solve_modelKARC(H, g, nlp_stop, X::PDataKARC, α::T) where {T}
    # target value should be close to satisfy αλ=||d||
    # target = [(abs(α * X.shifts[i] - X.norm_dirs[i])) for i = 1:X.nshifts]
    # positives = findall(X.positives)
    start = findfirst(X.positives)
    if isnothing(start)
        start = 1
    end
    positives = start:length(X.positives)
    target = ((abs(α * X.shifts[i] - X.norm_dirs[i])) for i in positives)

    # pick the closest shift to the target within positive definite H+λI
    indmin = argmin(target) # findmin(target[positives[start:end]])
    X.indmin = start + indmin - 1
    p_imin = X.indmin # Is there a diff with findall(X.positives)[X.indmin] ?
    X.d .= X.xShift[p_imin]
    X.λ = X.shifts[p_imin]

    return X.d, X.λ
end
