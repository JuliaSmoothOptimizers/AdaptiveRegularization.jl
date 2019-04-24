mutable struct TparamsST{T} <: Tparams{T}  # specific parameters for this Krylov variant
    τ :: T

    function TparamsST{T}() where T
        τ = 0.5
        return new{T}(τ)
    end
end
