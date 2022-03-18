function ARCqKOp(nlpstop::NLPStopping, shifts = 10.0 .^ (collect(-20.0:1.0:20.0)); kwargs...)
    return TRARC(
        nlpstop;
        c = Combi(
            HessOp,
            PDataKARC,
            solve_modelKARC,
        ),
        shifts = 10.0 .^ (collect(-20.0:1.0:20.0)),
        kwargs...,
    )
end
