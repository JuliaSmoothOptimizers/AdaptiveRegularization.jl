function ARCLDLt(nlpstop::NLPStopping; kwargs...)

    T = eltype(nlpstop.pb.meta.x0)

    return TRARC(
        nlpstop;
        TR = TrustRegion(T(10.0)),
        c = Combi(
            HessDense,
            PDataLDLt,
            solve_modelARCDiag,
            preprocessLDLt,
        ),
        kwargs...,
    )
end
