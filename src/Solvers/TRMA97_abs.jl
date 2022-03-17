function TRMA97_abs(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            HessDense,
            PDataMA97,
            solve_modelTRDiagAbs,
            preprocessMA97,
            Tparam{T}(),
        ),
        kwargs...,
    )
end
