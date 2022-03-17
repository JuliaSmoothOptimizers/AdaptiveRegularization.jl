function ARCMA57_abs(nlpstop::NLPStopping; kwargs...)
    return TRARC(
        nlpstop;
        c = Combi(
            HessSparse,
            PDataMA57,
            solve_modelARCDiagAbs,
        ),
        kwargs...,
    )
end
