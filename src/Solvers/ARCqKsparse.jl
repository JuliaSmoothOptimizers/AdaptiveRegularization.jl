function ARCqKsparse(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    shifts = 10.0 .^ (collect(-20.0:1.0:20.0))

    return TRARC(
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            HessSparse,
            PDataKARC,
            solve_modelKARC,
            preprocessKARC,
            TparamsKARC{T}(shifts),
        ),
        kwargs...,
    )
end
