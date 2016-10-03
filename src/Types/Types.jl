abstract TPData  # Ancestor of all PreProcess data

abstract PDataFact <: TPData  # Variants using matricial factorization

abstract PDataIter <: TPData  # Variants using iterative (Krylov) solvers

abstract Tparams  # Ancestor of variant specific parameters

type Tparam <: Tparams
end


type PDataK <: PDataIter
    d :: Array{Float64,1} # (H+λI)\g ; on first call = g
    λ :: Float64  # "active" value of λ; on first call = 0
    τ :: Float64  # Inexact Newton order parameter: stop when ||∇q||<||g||^(1+τ)

    indmin :: Int # index of best shift value  within "positive". On first call = 0
    
    positives :: Array{Int,1}  # indices of the shift values yielding (H+λI)⪰0
    xShift :: Array{Float64,2} # solutions for each shifted system
    shifts :: Array{Float64,1} # values of the shifts
    nshifts :: Int# number of shifts
    norm_dirs :: Array{Float64,1} # norms of xShifts
    OK :: Bool    # error flag when false
end


type PDataST <: PDataIter
    H # Hessian representation, most likely as a linear operator or a sparse matrix
    g # gradient vector
    τ :: Float64  # Inexact Newton order parameter: stop when ||∇q||<||g||^(1+τ)
    
    OK :: Bool
end

type PDataLDLt <: PDataFact
    L::Array{Float64,2}  # could be sparse
    D::Array{Float64,2}  # block diagonal 1X1 and 2X2
    P::Array{Float64,2}  # permutation matrix: should be vector
                         # LDL' = PHP'
    Δ::Array{Float64,1}  # diagonal, eigenvalues of D
    Q::Array{Float64,2}  # orthogonal matrix, eigenvectors of D:  should be sparse
                         # QΔQ'  =  D
    g̃::Array{Float64,1}  # transformed gradient 
    λ::Float64
    success::Bool        # previous iteration was successfull
    OK::Bool

    PDataLDLt() = new()
    PDataLDLt(L,D,P,DiagD,Q,g,l,success,OK) = new(L,D,P,DiagD,Q,g,l,success,OK)
end

type PDataSpectral <: PDataFact
    V::Array{Float64,2}   # orthogonal matrix, eigenvectors of Hessian matrix
    Δ::Array{Float64,1}   # eigenvalues of Hessian matrix
    g̃::Array{Float64,1}   # transformed gradient 
    λ::Float64
    success::Bool
    OK::Bool

    PDataSpectral() = new()
    PDataSpectral(V,Λ,g,l,success,OK) = new(V,Λ,g,l,success,OK)
end
