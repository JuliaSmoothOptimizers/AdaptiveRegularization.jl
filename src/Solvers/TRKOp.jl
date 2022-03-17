function TRKOp(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        c = Combi(
            HessOp,
            PDataKTR,
            solve_modelKTR,
        ),
        kwargs...,
    )
end
