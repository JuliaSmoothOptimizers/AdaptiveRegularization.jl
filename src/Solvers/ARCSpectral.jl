function ARCSpectral(nlpstop::NLPStopping; kwargs...)
    return TRARC(
        nlpstop;
        c = Combi(
            ARCTR.HessDense,
            ARCTR.PDataSpectral,
            ARCTR.solve_modelARCDiag,
        ),
        kwargs...,
    )
end
