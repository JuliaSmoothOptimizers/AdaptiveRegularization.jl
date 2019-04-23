export preprocessLDLt, preprocessLDLt_Shamanskii
function preprocessLDLt(H ,g, params :: Tparams, n1, n2)
    n, = size(g)
    global L = Matrix{Float64}(undef, n, n)
    global D = Matrix{Float64}(undef, n, n)
    global pp = Vector{Int}(undef, n)
    global ρ = nothing
    global ncomp = nothing

    try
        (L, D, pp, ρ, ncomp) = ldlt_symm(H, 'r')
    catch
 	println("*******   Problem in LDLt")
        res = PDataLDLt()
        res.OK = false
        return res
    end

    if true in isnan.(D)
 	println("*******   Problem in D from LDLt: NaN")
        println(" cond (H) = $(cond(H))")
        res = PDataLDLt()
        res.OK = false
        return res
    end

    X = eigen(D)
    Δ = X.values
    Q =  X.vectors
    l_m, = findmin(Δ)
    ĝ = L \ (g[pp])
    g̃ = Q' * ĝ
    n_g = norm(g)
    λ =  max(-l_m, 0.0)
    return  PDataLDLt(L, D, pp, Δ, Q, g̃, λ, true, true)
end


function AInv(X :: PDataLDLt, d̃ ::  Array{Float64,1})
    d̂ = X.Q * d̃
    u = X.L' \ d̂
    return u[invperm(X.pp)]
end


function reconstructH(X :: PDataLDLt)
    A = X.L * X.D * X.L'
    return A[X.pp.X.pp]
end
