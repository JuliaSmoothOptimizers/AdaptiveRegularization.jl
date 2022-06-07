export PDataKTR, PDataKARC, PDataST, PDataSpectral

abstract type TPData{T} end  # Ancestor of all PreProcess data

abstract type PDataFact{T} <: TPData{T} end # Variants using matricial factorization

abstract type PDataIter{T} <: TPData{T} end # Variants using iterative (Krylov) solvers

function preprocess(PData::TPData, H, g, n1, n2, α)
    return PData
end

mutable struct PDataKARC{T} <: PDataIter{T}
    d::Array{T,1}             # (H+λI)\g ; on first call = g
    λ::T                      # "active" value of λ; on first call = 0
    ζ::T                      # Inexact Newton order parameter: stop when ||∇q|| < ξ * ||g||^(1+ζ)
    ξ::T                      # Inexact Newton order parameter: stop when ||∇q|| < ξ * ||g||^(1+ζ)
    maxtol::T                 # Largest tolerance for Inexact Newton
    mintol::T                 # Smallest tolerance for Inexact Newton
    cgatol::Any
    cgrtol::Any

    indmin::Int               # index of best shift value  within "positive". On first call = 0

    positives::Array{Bool,1}   # indices of the shift values yielding (H+λI)⪰0
    xShift::Array{Array{T,1},1}        # solutions for each shifted system
    shifts::Array{T,1}        # values of the shifts
    nshifts::Int              # number of shifts
    norm_dirs::Array{T,1}     # norms of xShifts
    OK::Bool                  # preprocess success

    solver::CgLanczosShiftSolver
end

function PDataKARC(
    ::Type{S},
    ::Type{T},
    n;
    ζ = T(0.5),
    ξ = T(0.01),
    maxtol = T(0.01),
    mintol = sqrt(eps(T)),
    cgatol = (ζ, ξ, maxtol, mintol, gNorm2) -> max(mintol, min(maxtol, ξ * gNorm2^(1 + ζ))),
    cgrtol = (ζ, ξ, maxtol, mintol, gNorm2) -> max(mintol, min(maxtol, ξ * gNorm2^ζ)),
    shifts = 10.0 .^ collect(-20.0:1.0:20.0),
    kwargs...,
) where {S,T}
    d = S(undef, n)
    λ = zero(T)
    indmin = 1
    nshifts = length(shifts)
    positives = Array{Bool,1}(undef, nshifts)
    xShift = Array{S,1}(undef, nshifts)
    for i = 1:nshifts
        xShift[i] = S(undef, n)
    end
    norm_dirs = S(undef, nshifts)
    OK = true
    solver = CgLanczosShiftSolver(n, n, nshifts, S)
    return PDataKARC(
        d,
        λ,
        ζ,
        ξ,
        maxtol,
        mintol,
        cgatol,
        cgrtol,
        indmin,
        positives,
        xShift,
        T.(shifts),
        nshifts,
        norm_dirs,
        OK,
        solver,
    )
end

mutable struct PDataKTR{T} <: PDataIter{T}
    d::Array{T,1}             # (H+λI)\g ; on first call = g
    λ::T                      # "active" value of λ; on first call = 0
    ζ::T                      # Inexact Newton order parameter: stop when ||∇q|| < ξ * ||g||^(1+ζ)
    ξ::T                      # Inexact Newton order parameter: stop when ||∇q|| < ξ * ||g||^(1+ζ)
    maxtol::T                 # Largest tolerance for Inexact Newton
    mintol::T                 # Smallest tolerance for Inexact Newton
    cgatol::Any
    cgrtol::Any

    indmin::Int               # index of best shift value  within "positive". On first call = 0

    positives::Array{Bool,1}   # indices of the shift values yielding (H+λI)⪰0
    xShift::Array{Array{T,1},1}        # solutions for each shifted system
    shifts::Array{T,1}        # values of the shifts
    nshifts::Int              # number of shifts
    norm_dirs::Array{T,1}     # norms of xShifts
    OK::Bool                  # preprocess success

    solver::CgLanczosShiftSolver
end

function PDataKTR(
    ::Type{S},
    ::Type{T},
    n;
    ζ = T(0.5),
    ξ = T(0.01),
    maxtol = T(0.01),
    mintol = sqrt(eps(T)),
    cgatol = (ζ, ξ, maxtol, mintol, gNorm2) -> max(mintol, min(maxtol, ξ * gNorm2^(1 + ζ))),
    cgrtol = (ζ, ξ, maxtol, mintol, gNorm2) -> max(mintol, min(maxtol, ξ * gNorm2^ζ)),
    shifts = T[0.0; 10.0 .^ (collect(-20.0:1.0:20.0))],
    kwargs...,
) where {S,T}
    d = S(undef, n)
    λ = zero(T)
    indmin = 1
    nshifts = length(shifts)
    positives = Array{Bool,1}(undef, nshifts)
    xShift = Array{S,1}(undef, nshifts)
    for i = 1:nshifts
        xShift[i] = S(undef, n)
    end
    norm_dirs = S(undef, nshifts)
    OK = true
    solver = CgLanczosShiftSolver(n, n, nshifts, S)
    return PDataKTR(
        d,
        λ,
        ζ,
        ξ,
        maxtol,
        mintol,
        cgatol,
        cgrtol,
        indmin,
        positives,
        xShift,
        shifts,
        nshifts,
        norm_dirs,
        OK,
        solver,
    )
end

mutable struct PDataST{S,T} <: PDataIter{T}
    d::S
    λ::T
    ζ::T                      # Inexact Newton order parameter: stop when ||∇q|| < ξ * ||g||^(1+ζ)
    ξ::T                      # Inexact Newton order parameter: stop when ||∇q|| < ξ * ||g||^(1+ζ)
    maxtol::T                 # Largest tolerance for Inexact Newton
    mintol::T                 # Smallest tolerance for Inexact Newton
    cgatol::Any
    cgrtol::Any

    OK::Bool    # preprocess success
    solver::CgSolver
end

function PDataST(
    ::Type{S},
    ::Type{T},
    n;
    ζ = T(0.5),
    ξ = T(0.01),
    maxtol = T(0.01),
    mintol = sqrt(eps(T)),
    cgatol = (ζ, ξ, maxtol, mintol, gNorm2) -> max(mintol, min(maxtol, ξ * gNorm2^(1 + ζ))),
    cgrtol = (ζ, ξ, maxtol, mintol, gNorm2) -> max(mintol, min(maxtol, ξ * gNorm2^ζ)),
    kwargs...,
) where {S,T}
    d = S(undef, n)
    λ = zero(T)
    OK = true
    solver = CgSolver(n, n, S)
    return PDataST(d, λ, ζ, ξ, maxtol, mintol, cgatol, cgrtol, OK, solver)
end
