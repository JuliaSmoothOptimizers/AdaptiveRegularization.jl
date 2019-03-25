function preprocessLDLt(H ,g, params::Tparams,n1,n2)
    global L = Array(Float64,2)
    global D = Array(Float64,2)
    global pp = Array(Int,1)
    global ρ = Float64
    global ncomp = Int64

    try
        (L, D, pp, rho, ncomp) = ldlt_symm(H,'r')
    catch
 	println("*******   Problem in LDLt")
        res = PDataLDLt()
        res.OK = false
        return res
    end

    if true in isnan(D)
 	println("*******   Problem in D from LDLt: NaN")
        println(" cond (H) = $(cond(H))")
        res = PDataLDLt()
        res.OK = false
        return res
    end

    Δ, Q = eig(D)
    l_m, = findmin(Δ)
    ĝ = L\(g[pp])
    g̃ = Q'*ĝ
    n_g = norm(g)
    λ =  max(-l_m,0.0)
    return  PDataLDLt(L,D,pp,Δ,Q,g̃,λ,true,true)
end


function AInv(X :: PDataLDLt, d̃ ::  Array{Float64,1})
    d̂ = X.Q*d̃
    u = X.L'\d̂
    return u[invperm(X.pp)]
end


function reconstructH(X :: PDataLDLt)
    A = X.L*X.D*X.L'
    return A[X.pp.X.pp]
end
