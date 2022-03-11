mutable struct PDataMA97{T} <: PDataFact{T}
    PLD::Ma97          # H = (PL)D(PL)'  or  LDL' = P'HP
    Δ::Array{T,1}  # diagonal, eigenvalues of D
    D::Any
    H::Any
    pivot::Any
    #Q::Array{Float64,2}  # orthogonal matrix, eigenvectors of D:  should be sparse
    # QΔQ'  =  D
    Q::SparseMatrixCSC{T,Int64}
    g̃::Array{T,1}  # transformed gradient
    g::Array{T,1}  # untransformed gradient
    λ::T
    success::Bool        # previous iteration was successfull
    OK::Bool             # preprocess success

    PDataMA97() = new{Nothing}()
    PDataMA97(PLD, Δ, D, H, pivot, Q, g̃, g, l, success, OK) =
        new{eltype(g̃)}(PLD, Δ, D, H, pivot, Q, g̃, g, l, success, OK)
end


function preprocessMA97(H, g, params::Tparams, n1, n2)
    PLD = Ma97
    pivot = Array{Int32,1}
    vD1 = Array{Float64,2}
    vD2 = Array{Float64,2}

    H97 = convert(SparseMatrixCSC{Float64,Int}, H)
    try
        PLD = Ma97(H97, print_level = -1)
        ma97_factorize!(PLD)
    catch
        println("*******   Problem in MA97")
        res = PDataMA97()
        res.OK = false
        return res
    end
    try
        ret = ma97_inquire(PLD)
        vD = ret[2]
        pivot = ret[1]
        # Take care of zero valued pivots
        z(a) = a == 0
        pivot[findall(z, pivot)] = collect(1:length(pivot))[findall(z, pivot)]

        vD1 = copy(vD[1:1, :])'[:, 1:1]
        vD2 = copy(vD[2:2, 1:end-1])'[:, 1:1]
    catch
        println("*******   Problem after MA97")
        println(" Cond(H) = $(cond(full(H)))")
        res = PDataMA97()
        res.OK = false
        return res
    end

    vD1 = vec(vD1)
    vD2 = vec(vD2)
    D = SparseArrays.spdiagm(0 => vD1, -1 => vD2, 1 => vD2)

    #################  Future object BlockDiag operator?
    vQ1 = ones(length(vD1))       # vector representation of orthogonal matrix Q
    vQ2 = zeros(length(vD2))      #
    vQ2 = zeros(length(vD2))      #
    vQ2m = zeros(length(vD2))     #
    veig = copy(vD1)      # vector of eigenvalues of D, initialized to diagonal of D
    # if D diagonal, nothing more will be computed

    i = 1
    while i < length(vD1)
        if vD2[i] == 0.0
            i += 1
        else
            mA = [vD1[i] vD2[i]; vD2[i] vD1[i+1]]  #  2X2 submatrix
            # DiagmA, Qma = eig(mA)                 #  spectral decomposition of mA
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

    Q = sparse(SparseArrays.spdiagm(0 => vQ1, -1 => vQ2m, 1 => vQ2)) # sparse representation of Q

    Q = Q[abs.(pivot), abs.(pivot)]
    veig = veig[abs.(pivot)]

    Δ = 1.0 ./ veig #'[:,1]             # don't forget inquire returns D^-1
    l_m, = findmin(Δ)

    ĝ = copy(g)
    ma97_solve!(PLD, ĝ, job = :PL)

    g̃ = Q' * ĝ

    λ = max(-l_m, 0.0)
    return PDataMA97(PLD, Δ, D, H, pivot, Q, g̃, g, λ, true, true)
end

function AInv(X::PDataMA97, d̃::Array{T,1}) where {T}
    d̂ = X.Q * d̃

    ma97_solve!(X.PLD, d̂, job = :LPS)

    return d̂
end


# remove
function reconstructH(X::PDataMA97)
    return X.P' * X.L * X.D * X.L' * X.P
end
