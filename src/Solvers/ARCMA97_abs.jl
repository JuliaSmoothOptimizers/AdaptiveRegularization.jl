function ARCMA97_abs(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlp,
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            hessian_dense,
            PDataMA97{T},
            solve_modelARCDiagAbs,
            preprocessMA97,
            decreaseFact,
            Tparam{T}(),
        ),
        kwargs...,
    )
end
