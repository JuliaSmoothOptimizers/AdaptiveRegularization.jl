function ST_TROp(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            hessian_operator,
            PDataST{T},
            solve_modelST_TR,
            preprocessST,
            decreaseGen,
            TparamsST{T}(),
        ),
        kwargs...,
    )
end
