mutable struct TparamsST <: Tparams  # specific parameters for this Krylov variant
    τ :: Float64

    function TparamsST()
        τ = 0.5
        return new(τ)
    end
end
