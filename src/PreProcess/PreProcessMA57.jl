export PDataMA57, AInv

mutable struct PDataMA57{T} <: PDataFact{T}
    L::SparseMatrixCSC{T,Int64} # sparse L
    D::SparseMatrixCSC{T,Int64} # block diagonal 1X1 and 2X2
    #P::Array{T,2}  # permutation matrix: should be vector
    # LDL' = PHP'
    pp::Array{Int64,1}   # permutation vector LDL' = H[pp,pp]
    # P*v = v[pp]
    # P'*v = v[invperm(pp)]
    s::Array{T,1}  # Scaling matrix S=diagm(s)
    Δ::Array{T,1}  # diagonal, eigenvalues of D
    Q::SparseMatrixCSC{T,Int64}  # orthogonal matrix, eigenvectors of D:  should be sparse
    # QΔQ'  =  D
    g̃::Array{T,1}  # transformed gradient
    λ::T
    success::Bool        # previous iteration was successfull
    OK::Bool             # preprocess success

    PDataMA57() = new{Nothing}()
    PDataMA57(L, D, pp, s, Δ, Q, g̃, l, success, OK) =
        new{eltype(g̃)}(L, D, pp, s, Δ, Q, g̃, l, success, OK)
end

function preprocess(PData::PDataMA57, H, g, n1, n2)
    M = Ma57
    L = SparseMatrixCSC{Float64,Int64}
    D57 = SparseMatrixCSC{Float64,Int64}
    pp = Array{Int64,1}
    s = Array{Float64}
    ρ = Float64
    ncomp = Int64

    H57 = convert(SparseMatrixCSC{Float64,Int64}, H)  #  Hard coded Int
    try
        M = Ma57(H57)#, print_level = -1)
        ma57_factorize(M)
    catch
        println("*******   Problem in MA57_0")
        M = Ma57(H57, print_level = -1)
        ma57_factorize(M)
        res = PDataMA57()
        res.OK = false
        return res
    end
    ###### Tried in an alternative version ####
    # M.control.cntl[1] = 1.0 # diff
    # M.control.icntl[7] = 0 # diff
    ######################
    try
        (L, D57, s, pp) = ma57_get_factors(M)
    catch
        println("*******   Problem after MA57_1")
        # println(" Cond(H) = $(cond(full(H)))")
        res = PDataMA57_0()
        res.OK = false
        return res
    end

    #################  Future object BlockDiag operator?
    vD1 = diag(D57)       # create internal representation for block diagonal D
    vD2 = diag(D57, 1)     #
    vQ1 = ones(length(vD1))       # vector representation of orthogonal matrix Q
    vQ2 = zeros(length(vD2))      #
    vQ2m = zeros(length(vD2))     #
    veig = copy(vD1)      # vector of eigenvalues of D, initialized to diagonal of D
    # if D diagonal, nothing more will be computed

    i = 1
    while i < length(vD1)
        if vD2[i] == 0.0
            i += 1
        else
            mA = [vD1[i] vD2[i]; vD2[i] vD1[i+1]] #  2X2 submatrix
            # DiagmA, Qma = eig(mA)                   #  spectral decomposition of mA
            X = eigen(mA)
            DiagmA = X.values
            Qma = X.vectors
            veig[i] = DiagmA[1]
            vQ1[i] = Qma[1, 1]
            vQ2[i] = Qma[1, 2]
            vQ2m[i] = Qma[2, 1]
            vQ1[i+1] = Qma[2, 2]
            veig[i+1] = DiagmA[2]
            i += 2
        end
    end

    Q = sparse(SparseArrays.spdiagm(0 => vQ1, 1 => vQ2m, -1 => vQ2))
    Δ = veig

    l_m, = findmin(Δ)
    sg = s .* g
    L = SparseMatrixCSC{Float64,Int64}(L)  #### very important, the \ command doesn't work with Int32
    ĝ = L \ (sg[pp])
    g̃ = Q' * ĝ

    n_g = norm(g)
    λ = max(-l_m, 0.0)
    return PDataMA57(L, D57, pp, s, Δ, Q, g̃, λ, true, true)
end

function AInv(X::PDataMA57, d̃::Array{T,1}) where {T}
    d̂ = X.Q * d̃
    u = X.L' \ d̂
    return u[invperm(X.pp)] .* X.s
end
