export preprocessLDLt
function preprocessLDLt(H ,g, params :: Tparams, n1, n2)
    T = eltype(H)
    n, = size(g)
    global L = Matrix{T}(undef, n, n)
    global D = Matrix{T}(undef, n, n)
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

    # printstyled("dans preprocessLDLt apres ldlt_symm eltype(L) = $(eltype(L)) \n", color = :cyan)

    if true in isnan.(D)
 	println("*******   Problem in D from LDLt: NaN")
        println(" cond (H) = $(cond(H))")
        res = PDataLDLt()
        res.OK = false
        return res
    end

    # printstyled("dans preprocessLDLt \n", color = :yellow)
    # @show D[1]
    # @show D[2]
    # @show D[3]
    # @show D
    # println(" ")
    D = Hermitian(D)
    X = eigen(D)
    # printstyled("on a X \n", color = :yellow)
    Δ = X.values
    Q =  X.vectors
    l_m, = findmin(Δ)
    ĝ = L \ (g[pp])
    g̃ = Q' * ĝ

    n_g = norm(g)
    λ =  max(-l_m, 0.0)
    # printstyled("on sort de preprocessLDLt \n", color = :bold)
    return  PDataLDLt(L, D, pp, Δ, Q, g̃, λ, true, true)
end


function AInv(X :: PDataLDLt, d̃ ::  Array{T,1}) where T
    d̂ = X.Q * d̃
    u = X.L' \ d̂
    return u[invperm(X.pp)]
end


function reconstructH(X :: PDataLDLt)
    A = X.L * X.D * X.L'
    return A[X.pp.X.pp]
end
