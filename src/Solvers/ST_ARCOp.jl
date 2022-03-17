function ST_ARCOp(nlpstop::NLPStopping; kwargs...)
    return TRARC(
        nlpstop;
        c = Combi(
            HessOp,
            PDataST,
            solve_modelST_ARC,
        ),
        kwargs...,
    )
end
