function ST_ARCsparse(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        c = Combi(
            HessSparse,
            PDataST,
            solve_modelST_ARC,
        ),
        kwargs...,
    )
end
