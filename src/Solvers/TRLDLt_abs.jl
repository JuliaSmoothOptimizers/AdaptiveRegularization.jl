function TRLDLt_abs(nlpstop::NLPStopping; kwargs...)
    return TRARC(
        nlpstop;
        c = Combi(
            HessDense,
            PDataLDLt,
            solve_modelTRDiagAbs,
        ),
        kwargs...,
    )
end
