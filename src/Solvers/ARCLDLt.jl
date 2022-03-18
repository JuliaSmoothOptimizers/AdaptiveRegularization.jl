function ARCLDLt(nlpstop::NLPStopping; kwargs...)
    return TRARC(
        nlpstop;
        c = Combi(
            HessDense,
            PDataLDLt,
            solve_modelARCDiag,
        ),
        kwargs...,
    )
end
