abstract type Tparams{T} end

export Tparam

mutable struct Tparam{T} <: Tparams{T} end

export TparamsST

mutable struct TparamsST{T} <: Tparams{T}  # specific parameters for this Krylov variant
    ζ::T

    function TparamsST{T}() where {T}
        ζ = 0.5
        return new{T}(ζ)
    end
end

export TparamsKTR

mutable struct TparamsKTR{T} <: Tparams{T}  # specific parameters for this Krylov variant
    ζ::T
    shifts::Array{T,1}
    nshifts::Int
    function TparamsKTR{T}(; ζin::T = 0.5) where {T}
        ζ = ζin
        shifts = [0.0; 10.0 .^ (collect(-15.0:1.0:15.0))]
        nshifts = size(shifts, 1)
        return new{T}(ζ, shifts, nshifts)
    end
end

export TparamsKARC

mutable struct TparamsKARC{T} <: Tparams{T}  # specific parameters for this Krylov variant
    ζ::T
    τ::T
    shifts::Array{T,1}
    nshifts::Int
    function TparamsKARC{T}(shifts; τin::T = 1.0, ζin::T = 0.5) where {T}
        τ = τin  #  temporary testing
        ζ = ζin  #  inexact Newton
        nshifts = size(shifts, 1)
        return new{T}(ζ, τ, shifts, nshifts)
    end
end
