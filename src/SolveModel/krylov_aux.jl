
"""Given a regularization parameter `regulα`, a vector `x` and a direction `d`,
   return `σ` > 0 such that

    h(σ) ~ argmin_{α>0} h(α) = q(x + α d) +  ‖x + αd‖³/3regulα

in the Euclidean norm.
"""
function to_minimum(A :: LinearOperator, b :: Array{T,1},
                               x :: Vector{T},
                               p :: Vector{T}, Ap :: Vector{T}, pAp :: T,
                               α :: T,
                               regulα :: T) where T
    regulα > 0 || error("regulα must be positive")
    #println("regulα = ",regulα)

    Ax = copy(A * x)
    xAx = dot(x, Ax)
    xAp = dot(x, Ap)
    bx = dot(b, x)
    bp = dot(b, p)
    h = α -> 0.5 * (xAx + 2.0 * xAp * α + α^2 * pAp) - bx - α * bp + norm(x + α * p)^3/(3 * regulα)
    h0 = h(0)

    #dh = α -> xAp + α*pAp - bp + norm(x+α*p)*(x+α*p)'*p
    xx = dot(x, x)
    xp = dot(x, p)
    pp = dot(p, p)
    dh = α -> xAp + α * pAp - bp + sqrt(xx + 2 * α * xp + α^2 * pp) * (xp + α * pp) / regulα


#println("h0 = ",h(0),"  dh0 = ",dh(0))
    # search for a positive dh
    α = regulα
    while dh(α)<0
        α *= 5
    end
#println("h(α) = ",h(α),"  dh(α) = ",dh(α), "  α = ",α)

    a = T(0.0)
    σ = T(α/2)
    @assert ((dh(a) < 0) & (dh(α) > 0))
    while (((α - a) / max(α, a)) > 1e-8) & (abs(dh(σ)) > 1e-8)
        if dh(σ) < 0
            a = σ
        else
            α = σ
        end
        σ = (a + α) / 2
#println("  dh(σ) = ",dh(σ),"  dh(α) = ",dh(α),"  dh(a) = ",dh(a), "  σ = ",σ, "  α = ",α, "  a = ",a)
    end
#println("h(σ) = ",h(σ),"  dh(σ) = ",dh(σ), "  σ = ",σ)

    return σ
end
