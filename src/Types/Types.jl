abstract type TPData end  # Ancestor of all PreProcess data

abstract type PDataFact <: TPData  end # Variants using matricial factorization

abstract type PDataIter <: TPData  end # Variants using iterative (Krylov) solvers

abstract type Tparams  end # Ancestor of variant specific parameters

mutable struct Tparam <: Tparams end



mutable struct PDataK <: PDataIter
    d :: Array{Float64,1}      # (H+λI)\g ; on first call = g
    λ :: Float64               # "active" value of λ; on first call = 0
    τ :: Float64               # Inexact Newton order parameter: stop when ||∇q||<||g||^(1+τ )

    indmin :: Int              # index of best shift value  within "positive". On first call = 0

    positives :: Array{Int,1}  # indices of the shift values yielding (H+λI)⪰0
    xShift :: Array{Float64,2} # solutions for each shifted system
    shifts :: Array{Float64,1} # values of the shifts
    nshifts :: Int             # number of shifts
    norm_dirs :: Array{Float64,1} # norms of xShifts
    OK :: Bool                 # preprocess success
end


mutable struct PDataST <: PDataIter
    H             # Hessian representation, most likely as a linear operator or a sparse matrix
    g             # gradient vector
    τ :: Float64  # Inexact Newton order parameter: stop when ||∇q||<||g||^(1+τ)

    OK :: Bool    # preprocess success
end

mutable struct PDataLDLt <: PDataFact
    L        :: Array{Float64,2}   # could be sparse
    D        :: Array{Float64,2}   # block diagonal 1X1 and 2X2
    pp       :: Array{Int,1}       # permutation vector LDL' = H[pp,pp]
                                   # P*v = v[pp]
                                   # P'*v = v[invperm(pp)]
    #s::Array{Float64,1}           # Scaling matrix S=diagm(s)
    Δ        :: Array{Float64,1}   # diagonal, eigenvalues of D
    Q        :: Array{Float64,2}   # orthogonal matrix, eigenvectors of D:  should be sparse
                                   # QΔQ'  =  D
    g̃        ::Array{Float64,1}    # transformed gradient
    λ        ::Float64
    success  ::Bool                # previous iteration was successfull
    OK       ::Bool                # preprocess success

    PDataLDLt() = new()
    PDataLDLt(L, D, pp, Δ, Q, g, l, success, OK) = new(L, D, pp, Δ, Q, g, l, success, OK)
end

mutable struct PDataSpectral <: PDataFact
    V::Array{Float64,2}  # orthogonal matrix, eigenvectors of Hessian matrix
    Δ::Array{Float64,1}  # eigenvalues of Hessian matrix
    g̃::Array{Float64,1}  # transformed gradient
    λ::Float64
    success::Bool        # previous iteration was successfull
    OK::Bool             # preprocess success

    PDataSpectral() = new()
    PDataSpectral(V, Λ, g, l, success, OK) = new(V, Λ, g, l, success, OK)
end
