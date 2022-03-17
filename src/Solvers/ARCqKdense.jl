function ARCqKdense(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)

    return TRARC(
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            HessDense,
            PDataKARC,
            solve_modelKARC,
            preprocessKARC,
        ),
        shifts = 10.0 .^ (collect(-20.0:1.0:20.0)),
        kwargs...,
    )
end
