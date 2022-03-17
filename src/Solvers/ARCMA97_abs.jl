function ARCMA97_abs(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlp,
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            HessDense,
            PDataMA97{T},
            solve_modelARCDiagAbs,
            preprocessMA97,
            Tparam{T}(),
        ),
        kwargs...,
    )
end
