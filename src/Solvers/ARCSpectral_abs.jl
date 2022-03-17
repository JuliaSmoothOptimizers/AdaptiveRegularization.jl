function ARCSpectral_abs(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        c = Combi(
            ARCTR.HessDense,
            ARCTR.PDataSpectral,
            ARCTR.solve_modelARCDiagAbs,
        ),
        kwargs...,
    )
end
