export solve_modelKTR
function solve_modelKTR(nlp_stop, X :: PDataK, α:: T) where T
    # target value should be close to satisfy α=||d||
    target = [ ( abs( α - X.norm_dirs[i]) )   for i = 1 : X.nshifts ];

    # pick the closest shift to the target within positive definite H+λI
    # before, check the shift = 0 for direction within the trust region

    if (X.indmin ==0)  #  try Newton's direction
        X.indmin = 1
        if (X.positives[1]==1) & (X.norm_dirs[1] <= α)
            X.d = X.xShift[:,1]
            X.λ = 0.0
            return X.d, NaN * ones(length((X.d))), X.λ
        end
    end
    start = X.indmin
    bidon,indmin = findmin(target[X.positives[start:end]])
    X.indmin = start + indmin - 1
    p_imin = X.positives[X.indmin]
    X.d = X.xShift[:,p_imin]
    # Refine this rough constrained direction

    X.λ = X.shifts[p_imin]
    #println("In SolveModel   p_imin = $(p_imin),  α = $α")
    #println("In SolveModel   $(X.Ndirs[p_imin]-α) , $(X.Ndirs[p_imin+1]-α)  ,  α = $α")

    #end
    return X.d, NaN * ones(length((X.d))), X.λ
end

# To replace.
#    1. Form the signed ( α - X.norm_dirs[i]) ) within positives of the form +++...----
#    2. If starts with a +, interpolate the interval where sign changes
#       else (hard case) obtain from Krylov shift solver negative curvature directions
#            and extrapolate to the boundary the short (longest) direction
#            in the descent negative curvature direction
#
# Keep on treating λ=0 special for Newton's direction (no extrapolation to boundary)
