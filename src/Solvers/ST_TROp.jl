function ST_TROp(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            HessOp,
            PDataST,
            solve_modelST_TR,
            preprocessST,
        ),
        kwargs...,
    )
end
