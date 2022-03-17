function ARCMA57_Sham_λ(nlpstop::NLPStopping; λfact::Float64 = 1.0, kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            HessSparse,
            PDataMA57{T},
            (x, y, z) -> solve_modelARCDiag_HO(
                x,
                y,
                z,
                ho_correction = :Shamanskii_MA57,
                λfact = λfact,
            ),
            preprocessMA57,
            Tparam{T}(),
        ),
        kwargs...,
    )
end
