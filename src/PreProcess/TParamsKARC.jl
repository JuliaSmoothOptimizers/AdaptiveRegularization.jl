mutable struct TparamsKARC{T} <: Tparams{T}  # specific parameters for this Krylov variant
    ζ :: T
    τ :: T
    shifts :: Array{T, 1}
    nshifts :: Int
    function TparamsKARC{T}(shifts; τin::T = 1.0, ζin ::T = 0.5) where T
        τ = τin  #  temporary testing
        ζ = ζin  #  inexact Newton
        nshifts = size(shifts, 1)
        return new{T}(ζ, τ, shifts, nshifts)
    end
end
