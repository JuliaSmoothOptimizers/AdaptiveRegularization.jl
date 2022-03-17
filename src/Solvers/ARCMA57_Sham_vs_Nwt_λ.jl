function ARCMA57_Sham_vs_Nwt_λ(
    nlpstop::NLPStopping;
    λfact::Float64 = 1.0,
    nwt_res_fact = 0.25,
    kwargs...,
)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            HessSparse,
            PDataMA57,
            (x, y, z) -> solve_modelARCDiag_HO_vs_Nwt(
                x,
                y,
                z,
                λfact = λfact,
                nwt_res_fact = nwt_res_fact,
                ho_correction = :Shamanskii_MA57,
            ),
            preprocessMA57,
        ),
        kwargs...,
    )
end
