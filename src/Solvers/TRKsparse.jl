function TRKsparse(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        c = Combi(
            HessSparse,
            PDataKTR,
            solve_modelKTR,
        ),
        kwargs...,
    )
end
