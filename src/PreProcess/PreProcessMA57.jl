type PDataMA57 <: PDataFact
    L::SparseMatrixCSC{Float64,Int32} # sparse L
    D::SparseMatrixCSC{Float64,Int32} # block diagonal 1X1 and 2X2
    #P::Array{Float64,2}  # permutation matrix: should be vector
                         # LDL' = PHP'
    pp::Array{Int32,1}   # permutation vector LDL' = H[pp,pp]
                         # P*v = v[pp]
                         # P'*v = v[invperm(pp)]
    s::Array{Float64,1}  # Scaling matrix S=diagm(s)
    Δ::Array{Float64,1}  # diagonal, eigenvalues of D
    Q::SparseMatrixCSC{Float64,Int64}  # orthogonal matrix, eigenvectors of D:  should be sparse
                         # QΔQ'  =  D
    g̃::Array{Float64,1}  # transformed gradient 
    λ::Float64
    success::Bool        # previous iteration was successfull
    OK::Bool             # preprocess success

    PDataMA57() = new()
    PDataMA57(L,D,pp,s,Δ,Q,g,l,success,OK) = new(L,D,pp,s,Δ,Q,g,l,success,OK)
end

function preprocessMA57(H ,g, params::Tparams,n1,n2) 
    M = Ma57
    L = SparseMatrixCSC{Float64,Int32}
    D57 = SparseMatrixCSC{Float64,Int32}
    pp = Array(Int32,1)
    s = Array{Float64}
    ρ = Float64
    ncomp = Int64
    
    H57 = convert(SparseMatrixCSC{Cdouble,Int32}, H)  #  Hard coded Cdouble
    try
        M = Ma57(H,print_level=-1)
        ma57_factorize(M)
    catch
 	println("*******   Problem in MA57_0")
        M = Ma57(H,print_level=-1)
        ma57_factorize(M)        
        res = PDataMA57()
        res.OK = false
        return res
    end

    try
        (L, D57, s, pp) = ma57_get_factors(M)
    catch
        println("*******   Problem after MA57_0")
        println(" Cond(H) = $(cond(full(H)))")
        res = PDataMA57_0()
        res.OK = false
        return res
    end

    #################  Future object BlockDiag operator?
    vD1 = diag(D57)       # create internal representation for block diagonal D
    vD2 = diag(D57,1)     #

    vQ1 = ones(vD1)       # vector representation of orthogonal matrix Q
    vQ2 = zeros(vD2)      #
    vQ2 = zeros(vD2)      #
    vQ2m = zeros(vD2)     #
    veig = copy(vD1)      # vector of eigenvalues of D, initialized to diagonal of D
                          # if D diagonal, nothing more will be computed
    
    i=1;
    while i<length(vD1)
        if vD2[i] == 0.0
            i += 1
        else
            mA = [vD1[i] vD2[i];vD2[i] vD1[i+1]]  #  2X2 submatrix
            DiagmA, Qma = eig(mA)                 #  spectral decomposition of mA
            veig[i] = DiagmA[1]
            vQ1[i] = Qma[1,1]
            vQ2[i] = Qma[1,2]
            vQ2m[i] = Qma[2,1]
            vQ1[i+1] = Qma[2,2]
            veig[i+1] = DiagmA[2]
            i += 2
        end  
    end

    Q = spdiagm((vQ1,vQ2m,vQ2),[0,-1,1])           # sparse representation of Q
    
    Δ = veig

    l_m, = findmin(Δ)
    sg = s .* g
    ĝ = L\(sg[pp]) 
    g̃ = Q'*ĝ


    n_g = norm(g)
    λ =  max(-l_m,0.0) 
    return  PDataMA57(L,D57,pp,s,Δ,Q,g̃,λ,true,true)
end


function AInv(X :: PDataMA57, d̃ ::  Array{Float64,1})
    d̂ = X.Q*d̃
    u = X.L'\d̂
    return u[invperm(X.pp)] .* X.s
end


function reconstructH(X :: PDataMA57)
    A = X.L*X.D*X.L'
    return A[X.pp.X.pp]
end
