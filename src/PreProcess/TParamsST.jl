export TparamsST

mutable struct TparamsST{T} <: Tparams{T}  # specific parameters for this Krylov variant
    ζ :: T

    function TparamsST{T}() where T
        ζ = 0.5
        return new{T}(ζ)
    end
end
