function ARCMA57_Sham_vs_Nwt(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            HessSparse,
            PDataMA57,
            solve_modelARCDiag,
            preprocessMA57,
            Tparam{T}(),
        ),
        kwargs...,
    )
end
