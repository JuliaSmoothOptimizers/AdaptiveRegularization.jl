function ARCMA57(nlpstop::NLPStopping; kwargs...)
    return TRARC(
        nlpstop;
        c = Combi(
            HessSparse,
            PDataMA57,
            solve_modelARCDiag,
        ),
        kwargs...,
    )
end
