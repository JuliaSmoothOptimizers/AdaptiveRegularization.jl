mutable struct TparamsKARC{T} <: Tparams{T}  # specific parameters for this Krylov variant
    τ :: T
    shifts :: Array{T, 1}
    nshifts :: Int
    function TparamsKARC{T}(shifts) where T
        τ = 0.5
        nshifts = size(shifts, 1)
        return new{T}(τ, shifts, nshifts)
    end
end
