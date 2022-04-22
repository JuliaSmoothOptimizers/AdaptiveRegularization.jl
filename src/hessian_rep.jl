struct HessDense
    function HessDense(::Type{T}, ::Type{S}, n) where {T,S}
        return new()
    end
end
struct HessSparse
    function HessSparse(::Type{T}, ::Type{S}, n) where {T,S}
        return new()
    end
end
struct HessOp{S}
    Hv::S
    function HessOp(::Type{T}, ::Type{S}, n) where {T,S}
        return new{S}(S(undef, n))
    end
end
export HessDense, HessSparse, HessOp

function hessian!(workspace::HessDense, nlp, x)
    H = Matrix(hess(nlp, x))
    return H
end

function hessian!(workspace::HessOp, nlp, x)
    return hess_op!(nlp, x, workspace.Hv)
end

function hessian!(workspace::HessSparse, nlp, x)
    H = hess(nlp, x)
    return H
end
