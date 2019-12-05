export TparamsKTR

mutable struct TparamsKTR{T} <: Tparams{T}  # specific parameters for this Krylov variant
    ζ :: T
    shifts :: Array{T, 1}
    nshifts :: Int
    function TparamsKTR{T}(; ζin ::T = 0.5) where T
        ζ = ζin
        shifts = [0.0; 10.0.^(collect(-15.0:1.0:15.0))]
        nshifts = size(shifts, 1)
        return new{T}(ζ, shifts, nshifts)
    end
end
