function ARCMA57_abs(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            HessSparse,
            PDataMA57,
            solve_modelARCDiagAbs,
            preprocessMA57,
        ),
        kwargs...,
    )
end
