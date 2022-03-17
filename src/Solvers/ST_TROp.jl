function ST_TROp(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
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
