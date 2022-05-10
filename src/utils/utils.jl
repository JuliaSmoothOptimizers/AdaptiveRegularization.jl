export TrustRegion, convert_TR, convert_TR!, extract

"Exception type raised in case of error."
mutable struct TrustRegionException <: Exception
    msg::String
end

mutable struct TrustRegion{T}
    α₀::T
    α::T
    max_α::T
    acceptance_threshold::T
    increase_threshold::T
    reduce_threshold::T
    increase_factor::T
    decrease_factor::T
    max_unsuccinarow::Int

    function TrustRegion(
        α₀::T;
        max_α::T = T(1.0) / sqrt(eps(T)),
        acceptance_threshold::T = T(0.1),
        increase_threshold::T = T(0.75),
        reduce_threshold::T = T(0.1),
        increase_factor::T = T(5.0),
        decrease_factor::T = T(0.1),
        max_unsuccinarow::Int = 30,
    ) where {T}

        α₀ > T(0.0) || (α₀ = T(1.0))
        max_α > α₀ || throw(TrustRegionException("Invalid α₀"))
        (T(0.0) < acceptance_threshold < increase_threshold < T(1.0)) ||
            throw(TrustRegionException("Invalid thresholds"))
        (T(0.0) < decrease_factor < T(1.0) < increase_factor) ||
            throw(TrustRegionException("Invalid decrease/increase factors"))

        return new{T}(
            α₀,
            α₀,
            max_α,
            acceptance_threshold,
            increase_threshold,
            reduce_threshold,
            increase_factor,
            decrease_factor,
            max_unsuccinarow,
        )
    end
end


"""Compute the actual vs. predicted reduction radio ∆f/Δm, where
Δf = f - f_trial is the actual reduction is an objective/merit/penalty function,
Δq = q - q_trial is the reduction predicted by the model q of f.
We assume that q is being minimized, and therefore that Δq > 0.
"""
function compute_r(nlp, f, Δf, Δq, slope, d, xnext, gnext, robust)
    # If ever the next gradient is computed for round off errors reason, signal it.
    # # printstyled("dans compute r eltype(f) = $(eltype(f)) \n", color = :yellow)
    # # @show eltype(f)
    T = eltype(f)
    good_grad = false
    if robust & ((Δq < 10000 * eps(T)) | (abs(Δf) < 10000 * eps(T) * abs(f)))
        # trap potential cancellation errors
        grad!(nlp, xnext, gnext)
        good_grad = true
        slope_next = dot(gnext, d)

        Δf = -(slope_next + slope) / 2.0
    end
    r = Δf / Δq
    if isnan(r)
        r = 0.0
    end

    return r, good_grad, gnext
end

stop_norm(x) = norm(x, 2)
