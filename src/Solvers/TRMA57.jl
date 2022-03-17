function TRMA57(nlpstop::NLPStopping; kwargs...)

    return TRARC(
        nlpstop;
        c = Combi(
            HessSparse,
            PDataMA57,
            solve_modelTRDiag,
        ),
        kwargs...,
    )
end
