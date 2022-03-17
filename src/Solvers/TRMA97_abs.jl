function TRMA97_abs(nlpstop::NLPStopping; kwargs...)
    return TRARC(
        nlpstop;
        c = Combi(
            HessDense,
            PDataMA97,
            solve_modelTRDiagAbs,
        ),
        kwargs...,
    )
end
