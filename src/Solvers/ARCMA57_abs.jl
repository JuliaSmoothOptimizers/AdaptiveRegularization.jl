function ARCMA57_abs(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        c = Combi(
            HessSparse,
            PDataMA57,
            solve_modelARCDiagAbs,
        ),
        kwargs...,
    )
end
