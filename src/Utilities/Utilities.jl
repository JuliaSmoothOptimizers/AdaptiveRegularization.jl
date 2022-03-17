export TrustRegion, Combi, decreaseFact, convert_TR, convert_TR!, extract

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
    #params :: Tparams

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
        #params :: Tparams = Void)

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
        )#, params)
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


# A few headers displays
function print_header(Title)

    @printf(" \n\n %-15s\n", Title)
    @printf(
        "\n%-16s  %5s  %9s  %7s %7s  %5s  %5s  %5s  %6s  %s   %s\n",
        "Name",
        "nvar",
        "f",
        "‖∇f‖∞",
        "‖∇f‖₂",
        "#obj",
        "#grad",
        " H  ",
        "#hprod",
        "status",
        "time"
    )
end


function display_header_problems()
    @printf(
        "%-15s  %5s  %9s  %7s %7s  %5s  %5s  %6s  %s\n",
        "Name",
        "nvar",
        "f",
        "‖∇f‖∞",
        "‖∇f‖₂",
        "#obj",
        "#grad",
        "#hprod",
        "status"
    )
end

function display_header_iterations()
    @printf(" Iter    f          ||g||       λ        α       status \n",)
end

# A few displays for verbose iterations
function display_failure(iter, fnext, λ, α)
    @printf("%4d  %10.3e             %7.1e %7.1e    unsuccessful\n", iter, fnext, λ, α)
end

function display_v_success(iter, f, norm_g, λ, α)
    @printf("%4d  %10.3e %9.2e   %7.1e %7.1e Very successful\n", iter, f, norm_g, λ, α)
end

function display_success(iter, f, norm_g, λ, α)
    @printf("%4d  %10.3e %9.2e   %7.1e %7.1e      successful\n", iter, f, norm_g, λ, α)
end

function print_stats(prob, dim, f, gNorm, gnorm2, calls, status, timt)
    @printf(
        "%-16s  %5d  %9.2e  %7.1e  %7.1e  %5d  %5d %5d %6d  %s  %8.3f\n",
        prob,
        dim,
        f,
        gNorm,
        gnorm2,
        calls[1],
        calls[2],
        calls[3],
        calls[4],
        status,
        timt
    )
end



# default increase and decrease functions.
function decreaseBase(α::T, TR::TrustRegion) where {T}
    return α * TR.decrease_factor
end

function increase(α::T, TR::TrustRegion) where {T}
    return min(α * TR.increase_factor, TR.max_α)
end


function decreaseGen(X::TPData, α::T, TR::TrustRegion) where {T}
    return decreaseBase(α, TR)
end

function decreaseFact(X::PDataFact, α::T, TR::TrustRegion) where {T}
    X.success = false
    return decreaseBase(α, TR)
end


function increase(X::TPData, α::T, TR::TrustRegion) where {T}
    return increase(α, TR)
end


stop_norm(x) = norm(x, Inf)

# Valid combinations
#
mutable struct Combi{T, Hess, PData}
    solve_model::Function
    pre_process::Function
    decrease::Function
    params::Union{Tparams{T},Tparams}
    function Combi(
        ::Type{Hess},
        ::Type{PData},
        solve_model::Function,
        pre_process::Function,
        decrease::Function,
        params::Union{Tparams{T},Tparams},
    ) where {T, Hess, PData}
        return new{T, Hess, PData}(solve_model, pre_process, decrease, params)
    end
end

function extract(c::Combi)
    return c.solve_model, c.pre_process, c.decrease, c.params
end

function convert_TR(T, TR_init::TrustRegion)
    max_α_T = T(TR_init.max_α)
    acceptance_threshold_T = T(TR_init.acceptance_threshold)
    increase_threshold_T = T(TR_init.increase_threshold)
    reduce_threshold_T = T(TR_init.reduce_threshold)
    increase_threshold_T = T(TR_init.increase_threshold)
    decrease_factor_T = T(TR_init.decrease_factor)

    TRT = TrustRegion(
        T(TR_init.α),
        max_α = max_α_T,
        acceptance_threshold = acceptance_threshold_T,
        increase_threshold = increase_threshold_T,
        decrease_factor = decrease_factor_T,
    )

    return TRT
end
