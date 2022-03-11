function TRMA57(nlpstop::NLPStopping; kwargs...)

    return TRARC(
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            hessian_sparse,
            PDataMA57{T},
            solve_modelTRDiag,
            preprocessMA57,
            decreaseFact,
            Tparam{T}(),
        ),
        kwargs...,
    )
end
