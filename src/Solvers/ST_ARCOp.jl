function ST_ARCOp(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            HessOp,
            PDataST{T},
            solve_modelST_ARC,
            preprocessST,
            TparamsST{T}(),
        ),
        kwargs...,
    )
end
