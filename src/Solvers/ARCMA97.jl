function ARCMA97(nlpstop::NLPStopping; kwargs...)
    return TRARC(
        nlpstop;
        c = Combi(
            HessDense,
            PDataMA97,
            solve_modelARCDiag,
        ),
        kwargs...,
    )
end
