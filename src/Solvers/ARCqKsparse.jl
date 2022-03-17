function ARCqKsparse(nlpstop::NLPStopping; kwargs...)
    return TRARC(
        nlpstop;
        c = Combi(
            HessSparse,
            PDataKARC,
            solve_modelKARC,
        ),
        shifts = 10.0 .^ (collect(-20.0:1.0:20.0)),
        kwargs...,
    )
end
