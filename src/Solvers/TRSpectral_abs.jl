function TRSpectral_abs(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            HessDense,
            PDataSpectral{T},
            solve_modelTRDiagAbs,
            preprocessSpectral,
            decreaseFact,
            Tparam{T}(),
        ),
        kwargs...,
    )
end
