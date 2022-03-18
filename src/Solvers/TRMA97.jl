function TRMA97(nlpstop::NLPStopping; kwargs...)
    return TRARC(
        nlpstop;
        c = Combi(
            HessSparse,
            PDataMA97,
            solve_modelTRDiag,
        ),
        kwargs...,
    )
end
