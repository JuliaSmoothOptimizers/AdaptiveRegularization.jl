function ARCqKdense(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    shifts = 10.0 .^ (collect(-20.0:1.0:20.0))

    return TRARC(
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            HessDense,
            PDataKARC{T},
            solve_modelKARC,
            preprocessKARC,
            TparamsKARC{T}(shifts),
        ),
        kwargs...,
    )
end
