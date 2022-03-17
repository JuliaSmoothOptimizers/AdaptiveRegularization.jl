function TRMA57_Sham_vs_Nwt_λ(
    nlpstop::NLPStopping;
    corr_ho::Symbol = :Shamanskii_MA57,
    nwt_res_fact = 0.25,
    λfact = 1.0,
    kwargs...,
)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            HessSparse,
            PDataMA57{T},
            (x, y, z) -> solve_modelTRDiag_HO_vs_Nwt_λ(
                x,
                y,
                z;
                ho_correction = corr_ho,
                nwt_res_fact = nwt_res_fact,
                λfact = λfact,
            ),
            preprocessMA57,
            decreaseFact,
            Tparam{T}(),
        ),
        kwargs...,
    )
end
