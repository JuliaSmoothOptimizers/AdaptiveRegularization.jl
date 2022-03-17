function TRMA97(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            HessSparse,
            PDataMA97,
            solve_modelTRDiag,
            preprocessMA97,
            Tparam{T}(),
        ),
        kwargs...,
    )
end
