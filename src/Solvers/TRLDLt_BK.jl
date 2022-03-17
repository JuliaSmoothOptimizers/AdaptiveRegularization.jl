function TRLDLt_BK(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            HessDense,
            PDataLDLt,
            solve_modelTRDiag,
            preprocessLDLt2,
            Tparam{T}(),
        ),
        kwargs...,
    )
end
