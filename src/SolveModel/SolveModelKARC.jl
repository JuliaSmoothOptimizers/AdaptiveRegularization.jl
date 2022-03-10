function solve_modelKARC(nlp_stop, X :: PDataK, α:: T) where T
    # target value should be close to satisfy αλ=||d||
    target = [ ( abs( α*X.shifts[i] - X.norm_dirs[i]^X.τ) )   for i = 1 : X.nshifts ];

    # pick the closest shift to the target within positive definite H+λI
    X.indmin = max(1,X.indmin)
    start = X.indmin
    bidon,indmin = findmin(target[X.positives[start:end]])
    X.indmin = start + indmin - 1
    p_imin = X.positives[X.indmin]
    X.d = X.xShift[:,p_imin]
    X.λ = X.shifts[p_imin]

    return X.d, X.λ
end
