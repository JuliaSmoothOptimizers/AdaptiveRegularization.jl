function ST_TRdense(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            HessDense,
            PDataST,
            solve_modelST_TR,
            preprocessST,
            TparamsST{T}(),
        ),
        kwargs...,
    )
end
