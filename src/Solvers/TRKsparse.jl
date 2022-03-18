function TRKsparse(nlpstop::NLPStopping; kwargs...)
    return TRARC(
        nlpstop;
        c = Combi(
            HessSparse,
            PDataKTR,
            solve_modelKTR,
        ),
        kwargs...,
    )
end
