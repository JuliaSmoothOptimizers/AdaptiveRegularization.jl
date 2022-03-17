function ARCLDLt_abs(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        c = Combi(
            HessDense,
            PDataLDLt,
            solve_modelARCDiagAbs,
        ),
        kwargs...,
    )
end
