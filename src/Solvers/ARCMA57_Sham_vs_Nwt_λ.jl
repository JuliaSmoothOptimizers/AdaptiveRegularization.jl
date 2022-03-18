function ARCMA57_Sham_vs_Nwt_λ(
    nlpstop::NLPStopping;
    λfact::Float64 = 1.0,
    nwt_res_fact = 0.25,
    kwargs...,
)
    return TRARC(
        nlpstop;
        c = Combi(
            HessSparse,
            PDataMA57,
            (H, g, x, y, z) -> solve_modelARCDiag_HO_vs_Nwt(H, g, 
                x,
                y,
                z,
                λfact = λfact,
                nwt_res_fact = nwt_res_fact,
                ho_correction = :Shamanskii_MA57,
            ),
        ),
        kwargs...,
    )
end
