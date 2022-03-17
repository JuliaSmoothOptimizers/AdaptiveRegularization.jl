export preprocessLDLt, preprocessLDLt2
function preprocessLDLt(H, g, params::Tparams, n1, n2)
    # printstyled("dans preprocessLDLt ⤈ \n", color = :red)
    T = eltype(H)
    # @show T
    # @show H[1, 1]
    # @show H[1, 2]
    # @show H[2, 1]
    # @show H[2, 2]
    # @show g[1]
    # @show g[2]
    n, = size(g)
    # @show n
    global L = Matrix{T}(undef, n, n)
    global D = Matrix{T}(undef, n, n)
    global pp = Vector{Int}(undef, n)
    global ρ = nothing
    global ncomp = nothing

    # @show H

    try
        (L, D, pp, ρ, ncomp) = ldlt_symm(H, 'r')
    catch
        println("*******   Problem in LDLt")
        res = PDataLDLt()
        res.OK = false
        return res
    end

    # printstyled("dans preprocessLDLt apres ldlt_symm eltype(L) = $(eltype(L)) \n", color = :cyan)
    # @show L
    # @show D
    # @show L[1, 1]
    # @show L[2, 1]
    # @show L[1, 2]
    # @show L[2, 2]
    # @show D[1, 1]
    # @show D[2, 1]
    # @show D[1, 2]
    # @show D[2, 2]
    # @show pp
    # @show ρ
    # @show ncomp

    if true in isnan.(D)
        println("*******   Problem in D from LDLt: NaN")
        println(" cond (H) = $(cond(H))")
        res = PDataLDLt()
        res.OK = false
        return res
    end

    # @show L
    # @show D
    # @show pp

    # D = Hermitian(D)
    # @show D
    X = eigen(SymTridiagonal(D))
    Δ = X.values
    # @show Δ[1]
    # @show Δ[2]
    Q = X.vectors
    # @show Q[1, 1]
    # @show Q[2, 1]
    # @show Q[1, 2]
    # @show Q[2, 2]
    l_m, = findmin(Δ)
    # @show l_m
    ĝ = L \ (g[pp])
    # @show ĝ[1]
    # @show ĝ[2]
    # @show eltype(ĝ)
    # @show ĝ
    g̃ = Q' * ĝ
    # @show g̃[1]
    # @show g̃[2]
    # @show eltype(g̃)
    # @show g̃

    n_g = norm(g)
    # @show n_g
    λ = max(-l_m, 0.0)
    # @show eltype(λ)
    # @show eltype(λ)
    # @show λ
    # printstyled("avant de sortir de PreProcessLDLt ↑ \n", color = :red)
    return PDataLDLt(L, D, pp, Δ, Q, g̃, λ, true, true)
end

function preprocessLDLt2(H, g, params::Tparams, n1, n2)
    # printstyled("dans preprocessLDLt \n", color = :red)
    T = eltype(H)
    # @show H
    # @show g
    n, = size(g)
    # @show n
    global L = Matrix{T}(undef, n, n)
    global D = Matrix{T}(undef, n, n)
    global pp = Vector{Int}(undef, n)
    global ρ = nothing
    global ncomp = nothing
    global S = nothing

    # @show H
    # printstyled("ici ↩ \n", color = :red)

    try
        # (L, D, pp, ρ, ncomp) = ldlt_symm(H, 'r')
        # @show H
        S = bunchkaufman(Symmetric(H, :L), true, check = false)
    catch
        println("*******   Problem in LDLt")
        res = PDataLDLt()
        res.OK = false
        return res
    end
    L = S.L
    D = S.D
    pp = S.p

    # printstyled("dans preprocessLDLt apres ldlt_symm eltype(L) = $(eltype(L)) \n", color = :cyan)
    # @show L
    # @show D
    # @show pp
    # @show ρ
    # @show ncomp

    if true in isnan.(D)
        println("*******   Problem in D from LDLt: NaN")
        println(" cond (H) = $(cond(H))")
        res = PDataLDLt()
        res.OK = false
        return res
    end

    # @show L
    # @show D
    # @show pp

    D = Hermitian(D)
    X = eigen(Matrix(D))
    Δ = X.values
    # @show Δ
    Q = X.vectors
    l_m, = findmin(Δ)
    # @show l_m
    ĝ = L \ (g[pp])
    # @show ĝ[1]
    # @show ĝ[2]
    # @show eltype(ĝ)
    g̃ = Q' * ĝ
    # @show g̃[1]
    # @show g̃[2]
    # @show eltype(g̃)

    n_g = norm(g)
    # @show n_g
    λ = max(-l_m, 0.0)
    # @show λ
    # printstyled("avant de sortir de PreProcessLDLt ↑ \n", color = :red)
    return PDataLDLt(L, D, pp, Δ, Q, g̃, λ, true, true)
end

function AInv(X::PDataLDLt, d̃::Array{T,1}) where {T}
    d̂ = X.Q * d̃
    u = X.L' \ d̂
    # @show u[invperm(X.pp)][1]
    # @show u[invperm(X.pp)][2]
    # printstyled("avant la sortie de AInv ⟰ \n", color = :red)
    return u[invperm(X.pp)]
end
