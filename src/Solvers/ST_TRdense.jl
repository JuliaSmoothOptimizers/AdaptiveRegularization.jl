function ST_TRdense(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        c = Combi(
            HessDense,
            PDataST,
            solve_modelST_TR,
        ),
        kwargs...,
    )
end
