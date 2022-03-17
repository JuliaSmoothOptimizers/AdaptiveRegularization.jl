function ARCMA97_abs(nlpstop::NLPStopping; kwargs...)
    return TRARC(
        nlpstop;
        c = Combi(
            HessDense,
            PDataMA97,
            solve_modelARCDiagAbs,
        ),
        kwargs...,
    )
end
