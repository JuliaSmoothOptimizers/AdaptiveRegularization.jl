function ARCLDLt_HO_vs_Nwt(
    nlpstop::NLPStopping;
    corr_ho::Symbol = :Shamanskii,
    nwt_res_fact = 0.8,
    λfact = 100.0,
    kwargs...,
)
    return TRARC(
        nlpstop;
        c = Combi(
            HessDense,
            PDataLDLt,
            (H, g, x, y, z) -> solve_modelARCDiag_HO_vs_Nwt(H, g, 
                x,
                y,
                z,
                λfact = λfact,
                nwt_res_fact = nwt_res_fact,
            ),
        ),
        kwargs...,
    )
end
