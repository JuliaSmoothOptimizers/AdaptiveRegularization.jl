function TRKdense(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        c = Combi(
            HessDense,
            PDataKTR,
            solve_modelKTR,
        ),
        kwargs...,
    )
end
