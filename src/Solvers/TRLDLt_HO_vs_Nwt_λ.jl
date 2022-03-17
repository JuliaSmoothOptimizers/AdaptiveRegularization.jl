function TRLDLt_HO_vs_Nwt_λ(
    nlpstop::NLPStopping;
    corr_ho::Symbol = :Shamanskii,
    nwt_res_fact = 0.25,
    λfact = 1.0,
    kwargs...,
)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        TR = TrustRegion(T(10.0)),
        c = Combi(
            HessDense,
            PDataLDLt{T},
            (x, y, z) -> solve_modelTRDiag_HO_vs_Nwt_λ(
                x,
                y,
                z;
                ho_correction = corr_ho,
                nwt_res_fact = nwt_res_fact,
                λfact = λfact,
            ),
            preprocessLDLt,
            Tparam{T}(),
        ),
        kwargs...,
    )
end
