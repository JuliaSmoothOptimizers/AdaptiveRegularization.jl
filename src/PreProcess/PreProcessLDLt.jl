function preprocessLDLt(H ,g, params::Tparams,n1,n2)
    L = Array(Float64,2)
    D = Array(Float64,2)
    P = Array(Float64,2)
    ρ = Float64
    ncomp = Int64
    
    try
        (L, D, P, rho, ncomp) = ldlt_symm(H,'r')
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

    DiagD, Q = eig(D)
    l_m, = findmin(DiagD)
    ĝ = L\(P*g)   #(P'*L)\g
    g̃ = Q'*ĝ
    n_g = norm(g)
    λ =  max(-l_m,0.0) 
    return  PDataLDLt(L,D,P,DiagD,Q,g̃,λ,true,true)
end


function TtildeInv(X :: PDataLDLt, d̃ ::  Array{Float64,1})
    d̂ = X.Q*d̃
    return X.P'*(X.L'\d̂)
end


function reconstructH(X :: PDataLDLt)
    return X.P'*X.L*X.D*X.L'*X.P
end
