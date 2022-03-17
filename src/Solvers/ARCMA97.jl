function ARCMA97(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlp,
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            HessDense,
            PDataMA97,
            solve_modelARCDiag,
            preprocessMA97,
            Tparam{T}(),
        ),
        kwargs...,
    )
end
