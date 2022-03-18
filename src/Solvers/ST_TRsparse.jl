function ST_TRsparse(nlpstop::NLPStopping; kwargs...)
    return TRARC(
        nlpstop;
        c = Combi(
            HessSparse,
            PDataST,
            solve_modelST_TR,
        ),
        kwargs...,
    )
end
