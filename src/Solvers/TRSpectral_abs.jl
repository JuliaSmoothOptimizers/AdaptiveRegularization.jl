function TRSpectral_abs(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
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
