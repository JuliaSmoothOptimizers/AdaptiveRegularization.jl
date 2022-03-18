function TRSpectral(nlpstop::NLPStopping; kwargs...)
    return TRARC(
        nlpstop;
        c = Combi(
            HessDense,
            PDataSpectral,
            solve_modelTRDiag,
        ),
        kwargs...,
    )
end
