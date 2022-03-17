function ST_TRdense(nlpstop::NLPStopping; kwargs...)
    return TRARC(
        nlpstop;
        c = Combi(
            HessDense,
            PDataST,
            solve_modelST_TR,
        ),
        kwargs...,
    )
end
