function TRKOp(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            HessOp,
            PDataKTR,
            solve_modelKTR,
            preprocessKTR,
        ),
        kwargs...,
    )
end
