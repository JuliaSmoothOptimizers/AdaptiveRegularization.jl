function TRKdense(nlpstop::NLPStopping; kwargs...)
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
