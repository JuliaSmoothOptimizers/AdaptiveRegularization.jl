function ARCMA57(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            HessSparse,
            PDataMA57{T},
            solve_modelARCDiag,
            preprocessMA57,
            decreaseFact,
            Tparam{T}(),
        ),
        kwargs...,
    )
end
