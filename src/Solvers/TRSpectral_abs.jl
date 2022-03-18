function TRSpectral_abs(nlpstop::NLPStopping; kwargs...)
    return TRARC(
        nlpstop;
        c = Combi(
            HessDense,
            PDataSpectral,
            solve_modelTRDiagAbs,
        ),
        kwargs...,
    )
end
