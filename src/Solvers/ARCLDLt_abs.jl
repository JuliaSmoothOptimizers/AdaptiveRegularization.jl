function ARCLDLt_abs(nlpstop::NLPStopping; kwargs...)
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
