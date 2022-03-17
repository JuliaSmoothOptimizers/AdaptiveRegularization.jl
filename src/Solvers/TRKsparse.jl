function TRKsparse(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            HessSparse,
            PDataKTR,
            solve_modelKTR,
            preprocessKTR,
            TparamsKTR{T}(),
        ),
        kwargs...,
    )
end
