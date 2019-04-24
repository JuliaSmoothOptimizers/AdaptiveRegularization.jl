mutable struct TparamsKARC{T} <: Tparams{T}  # specific parameters for this Krylov variant
    τ :: T
    shifts :: Array{T, 1}
    nshifts :: Int
    function TparamsKARC{T}() where T
        τ = 0.5
        shifts = 10.0.^(collect(-15.0:1.0:15.0))
        nshifts = size(shifts, 1)
        return new{T}(τ, shifts, nshifts)
    end
end
