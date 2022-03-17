function ARCLDLt_HO_vs_Nwt(
    nlpstop::NLPStopping;
    corr_ho::Symbol = :Shamanskii,
    nwt_res_fact = 0.8,
    λfact = 100.0,
    kwargs...,
)

    T = eltype(nlpstop.pb.meta.x0)

    return TRARC(
        nlpstop;
        TR = TrustRegion(T(10.0)),
        c = Combi(
            HessDense,
            PDataLDLt,
            (x, y, z) -> solve_modelARCDiag_HO_vs_Nwt(
                x,
                y,
                z,
                λfact = λfact,
                nwt_res_fact = nwt_res_fact,
            ),
            preprocessLDLt,
        ),
        kwargs...,
    )
end
