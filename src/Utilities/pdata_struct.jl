export PDataIter
export PDataLDLt, PDataKTR, PDataKARC, PDataST

abstract type TPData{T} end  # Ancestor of all PreProcess data

abstract type PDataFact{T} <: TPData{T} end # Variants using matricial factorization

abstract type PDataIter{T} <: TPData{T} end # Variants using iterative (Krylov) solvers

mutable struct PDataKARC{T} <: PDataIter{T}
    d::Array{T,1}             # (H+λI)\g ; on first call = g
    λ::T                      # "active" value of λ; on first call = 0
    ζ::T                      # Inexact Newton order parameter: stop when ||∇q||<||g||^(1+ζ )

    indmin::Int               # index of best shift value  within "positive". On first call = 0

    positives::Array{Bool,1}   # indices of the shift values yielding (H+λI)⪰0
    xShift::Array{Array{T,1},1}        # solutions for each shifted system
    shifts::Array{T,1}        # values of the shifts
    nshifts::Int              # number of shifts
    norm_dirs::Array{T,1}     # norms of xShifts
    OK::Bool                  # preprocess success

    solver::CgLanczosShiftSolver
end

function PDataKARC(::Type{S}, ::Type{T}, n; ζ = 0.5, shifts = 10.0 .^ collect(-20.0:1.0:20.0), kwargs...) where {S, T}
    d = S(undef, n)
    λ = zero(T)
    indmin = 1
    nshifts = length(shifts)
    positives = Array{Bool, 1}(undef, nshifts)
    xShift = Array{S, 1}(undef, nshifts)
    for i=1:nshifts
        xShift[i] = S(undef, n)
    end
    norm_dirs = S(undef, nshifts)
    OK = true
    solver = CgLanczosShiftSolver(n, n, nshifts, S)
    return PDataKARC(d, λ, ζ, indmin, positives, xShift, shifts, nshifts, norm_dirs, OK, solver)
end

mutable struct PDataKTR{T} <: PDataIter{T}
    d::Array{T,1}             # (H+λI)\g ; on first call = g
    λ::T                      # "active" value of λ; on first call = 0
    ζ::T                      # Inexact Newton order parameter: stop when ||∇q||<||g||^(1+ζ )

    indmin::Int               # index of best shift value  within "positive". On first call = 0

    positives::Array{Bool,1}   # indices of the shift values yielding (H+λI)⪰0
    xShift::Array{Array{T,1},1}        # solutions for each shifted system
    shifts::Array{T,1}        # values of the shifts
    nshifts::Int              # number of shifts
    norm_dirs::Array{T,1}     # norms of xShifts
    OK::Bool                  # preprocess success

    solver::CgLanczosShiftSolver
end

function PDataKTR(::Type{S}, ::Type{T}, n; ζ = 0.5, shifts = [0.0; 10.0 .^ (collect(-20.0:1.0:20.0))], kwargs...) where {S, T}
    d = S(undef, n)
    λ = zero(T)
    indmin = 1
    nshifts = length(shifts)
    positives = Array{Bool, 1}(undef, nshifts)
    xShift = Array{S, 1}(undef, nshifts)
    for i=1:nshifts
        xShift[i] = S(undef, n)
    end
    norm_dirs = S(undef, nshifts)
    OK = true
    solver = CgLanczosShiftSolver(n, n, nshifts, S)
    return PDataKTR(d, λ, ζ, indmin, positives, xShift, shifts, nshifts, norm_dirs, OK, solver)
end

mutable struct PDataST{T} <: PDataIter{T}
    H::Any             # Hessian representation, most likely as a linear operator or a sparse matrix
    g::Any             # gradient vector
    ζ::T        # Inexact Newton order parameter: stop when ||∇q||<||g||^(1+ζ)

    OK::Bool    # preprocess success
end

function PDataST(::Type{S}, ::Type{T}, n; ζ = 0.5, kwargs...) where {S, T}
    OK = true
    return PDataST(nothing, nothing, ζ, OK)
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

function PDataLDLt(::Type{S}, ::Type{T}, n; kwargs...) where {S, T}
    L = Array{T, 2}(undef, n, n) # ::Array{T,2}          # could be sparse
    D = Array{T, 2}(undef, n, n) #::Array{T,2}          # block diagonal 1X1 and 2X2
    pp = collect(1:n) #::Array{Int,1}        # permutation vector LDL' = H[pp,pp]
    Δ = Array{T, 1}(undef, n) #::Array{T,1}          # diagonal, eigenvalues of D
    Q = Array{T, 2}(undef, n, n) #::Array{T,2}          # orthogonal matrix, eigenvectors of D:  should be sparse
    g̃ = Array{T, 1}(undef, n) #::Array{T,1}           # transformed gradient
    λ = zero(T) # ::T
    success = true::Bool                 # previous iteration was successfull
    OK = true
    return PDataLDLt(L, D, pp, Δ, Q, g̃, λ, success, OK)
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

function PDataSpectral(::Type{S}, ::Type{T}, n; kwargs...) where {S, T}
    V = Array{T, 2}(undef, n, n) # ::Array{T,2}          # could be sparse
    Δ = Array{T, 1}(undef, n) #::Array{T,1}          # diagonal, eigenvalues of D
    g̃ = Array{T, 1}(undef, n) #::Array{T,1}           # transformed gradient
    λ = zero(T) # ::T
    success = true::Bool                 # previous iteration was successfull
    OK = true
    return PDataSpectral(V, Δ, g̃, λ, success, OK)
end
