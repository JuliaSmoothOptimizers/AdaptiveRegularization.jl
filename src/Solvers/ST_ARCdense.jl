function ST_ARCdense(nlpstop::NLPStopping; kwargs...)
    return TRARC(
        nlpstop;
        c = Combi(
            HessDense,
            PDataST,
            solve_modelST_ARC,
        ),
        kwargs...,
    )
end
