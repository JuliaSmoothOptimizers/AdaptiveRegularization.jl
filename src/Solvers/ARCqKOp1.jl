function ARCqKOp1(nlpstop::NLPStopping; ζ = 0.5, τ = 1.0, kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    shifts = 10.0 .^ (collect(-15.0:1.0:15.0))
    return TRARC(
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            hessian_operator,
            PDataK{T},
            solve_modelKARC,
            preprocessKARC,
            decreaseKARC,
            TparamsKARC{T}(shifts, ζin = ζ, τin = τ),
        ),
        kwargs...,
    )
end
