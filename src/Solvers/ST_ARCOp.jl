function ST_ARCOp(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
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
