export PDataIter
export PDataLDLt, PDataK, PDataST

abstract type TPData{T} end  # Ancestor of all PreProcess data

abstract type PDataFact{T} <: TPData{T} end # Variants using matricial factorization

abstract type PDataIter{T} <: TPData{T} end # Variants using iterative (Krylov) solvers

mutable struct PDataK{T} <: PDataIter{T}
    d::Array{T,1}             # (H+λI)\g ; on first call = g
    λ::T                      # "active" value of λ; on first call = 0
    ζ::T                      # Inexact Newton order parameter: stop when ||∇q||<||g||^(1+ζ )
    τ::T                      # temporary testing parameter for decreaseARCqK

    indmin::Int               # index of best shift value  within "positive". On first call = 0

    positives::Array{Int,1}   # indices of the shift values yielding (H+λI)⪰0
    xShift::Array{Array{T,1},1}        # solutions for each shifted system
    shifts::Array{T,1}        # values of the shifts
    nshifts::Int              # number of shifts
    norm_dirs::Array{T,1}     # norms of xShifts
    OK::Bool                  # preprocess success
end


mutable struct PDataST{T} <: PDataIter{T}
    H::Any             # Hessian representation, most likely as a linear operator or a sparse matrix
    g::Any             # gradient vector
    ζ::T        # Inexact Newton order parameter: stop when ||∇q||<||g||^(1+ζ)

    OK::Bool    # preprocess success
end

mutable struct PDataLDLt{T} <: PDataFact{T}
    L::Array{T,2}          # could be sparse
    D::Array{T,2}          # block diagonal 1X1 and 2X2
    pp::Array{Int,1}        # permutation vector LDL' = H[pp,pp]
    # P*v = v[pp]
    # P'*v = v[invperm(pp)]
    #s::Array{T,1}                  # Scaling matrix S=diagm(s)
    Δ::Array{T,1}          # diagonal, eigenvalues of D
    Q::Array{T,2}          # orthogonal matrix, eigenvectors of D:  should be sparse
    # QΔQ'  =  D
    g̃::Array{T,1}           # transformed gradient
    λ::T
    success::Bool                 # previous iteration was successfull
    OK::Bool                 # preprocess success

    PDataLDLt() = new{Nothing}()
    PDataLDLt(L, D, pp, Δ, Q, g, l, success, OK) =
        new{eltype(L)}(L, D, pp, Δ, Q, g, l, success, OK)
end

mutable struct PDataSpectral{T} <: PDataFact{T}
    V::Array{T,2}           # orthogonal matrix, eigenvectors of Hessian matrix
    Δ::Array{T,1}           # eigenvalues of Hessian matrix
    g̃::Array{T,1}           # transformed gradient
    λ::T
    success::Bool           # previous iteration was successfull
    OK::Bool                # preprocess success

    PDataSpectral() = new{T}()
    PDataSpectral(V, Λ, g, l, success, OK) = new{eltype(V)}(V, Λ, g, l, success, OK)
end
