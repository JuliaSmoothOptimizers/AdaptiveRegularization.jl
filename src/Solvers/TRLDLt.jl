function TRLDLt(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        TR = TrustRegion(T(10.0)),
        c = Combi(
            HessDense,
            PDataLDLt{T},
            solve_modelTRDiag,
            preprocessLDLt,
            decreaseFact,
            Tparam{T}(),
        ),
        kwargs...,
    )
end
