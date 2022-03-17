function TRLDLt_HO_MP(nlpstop::NLPStopping; corr_ho::Symbol = :Shamanskii, kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC2_MP(
        nlp,
        nlpstop;
        TR = TrustRegion(10.0),
        c = Combi(
            HessDense,
            PDataLDLt,
            (x, y, z) -> solve_modelTRDiag_HO(x, y, z, ho_correction = corr_ho, fact = 2.0),
            preprocessLDLt,
        ),
        kwargs...,
    )
end
