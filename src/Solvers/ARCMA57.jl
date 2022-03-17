function ARCMA57(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        c = Combi(
            HessSparse,
            PDataMA57,
            solve_modelARCDiag,
        ),
        kwargs...,
    )
end
