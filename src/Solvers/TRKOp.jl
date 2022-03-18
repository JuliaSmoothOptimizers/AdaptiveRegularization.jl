function TRKOp(nlpstop::NLPStopping; kwargs...)
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
