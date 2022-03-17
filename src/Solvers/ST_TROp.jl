function ST_TROp(nlpstop::NLPStopping; kwargs...)
    return TRARC(
        nlpstop;
        c = Combi(
            HessOp,
            PDataST,
            solve_modelST_TR,
        ),
        kwargs...,
    )
end
