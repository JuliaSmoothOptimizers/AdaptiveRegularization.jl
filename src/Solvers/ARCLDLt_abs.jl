function ARCLDLt_abs(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            HessDense,
            PDataLDLt,
            solve_modelARCDiagAbs,
            preprocessLDLt,
            Tparam{T}(),
        ),
        kwargs...,
    )
end
