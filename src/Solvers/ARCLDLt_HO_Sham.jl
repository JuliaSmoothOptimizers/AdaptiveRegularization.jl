function ARCLDLt_HO_Sham(nlpstop::NLPStopping; λfact = 100.0, kwargs...)
    return TRARC(
        nlpstop;
        c = Combi(
            HessDense,
            PDataLDLt,
            (H, g, x, y, z) ->
                solve_modelARCDiag_HO(H, g, x, y, z, ho_correction = :Shamanskii, λfact = λfact),
        ),
        kwargs...,
    )
end
