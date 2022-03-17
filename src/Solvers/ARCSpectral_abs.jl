function ARCSpectral_abs(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            ARCTR.HessDense,
            ARCTR.PDataSpectral,
            ARCTR.solve_modelARCDiagAbs,
            ARCTR.preprocessSpectral,
            ARCTR.Tparam{T}(),
        ),
        kwargs...,
    )
end
