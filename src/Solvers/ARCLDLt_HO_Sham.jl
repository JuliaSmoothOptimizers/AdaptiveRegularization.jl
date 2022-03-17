function ARCLDLt_HO_Sham(nlpstop::NLPStopping; λfact = 100.0, kwargs...)

    T = eltype(nlpstop.pb.meta.x0)

    return TRARC(
        nlpstop;
        c = Combi(
            HessDense,
            PDataLDLt,
            (x, y, z) ->
                solve_modelARCDiag_HO(x, y, z, ho_correction = :Shamanskii, λfact = λfact),
        ),
        kwargs...,
    )
end
