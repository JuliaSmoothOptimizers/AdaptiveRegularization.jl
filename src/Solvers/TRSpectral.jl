function TRSpectral(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            HessDense,
            PDataSpectral,
            solve_modelTRDiag,
            preprocessSpectral,
            Tparam{T}(),
        ),
        kwargs...,
    )
end
