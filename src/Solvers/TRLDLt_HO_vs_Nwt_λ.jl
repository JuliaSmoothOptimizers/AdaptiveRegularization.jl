function TRLDLt_HO_vs_Nwt_λ(
    nlpstop::NLPStopping;
    corr_ho::Symbol = :Shamanskii,
    nwt_res_fact = 0.25,
    λfact = 1.0,
    kwargs...,
)
    return TRARC(
        nlpstop;
        c = Combi(
            HessDense,
            PDataLDLt,
            (H, g, x, y, z) -> solve_modelTRDiag_HO_vs_Nwt_λ(H, g, 
                x,
                y,
                z;
                ho_correction = corr_ho,
                nwt_res_fact = nwt_res_fact,
                λfact = λfact,
            ),
        ),
        kwargs...,
    )
end
