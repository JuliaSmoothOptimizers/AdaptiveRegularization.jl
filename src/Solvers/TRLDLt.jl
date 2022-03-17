function TRLDLt(nlpstop::NLPStopping; kwargs...)
    return TRARC(
        nlpstop;
        c = Combi(
            HessDense,
            PDataLDLt,
            solve_modelTRDiag,
        ),
    )
end
