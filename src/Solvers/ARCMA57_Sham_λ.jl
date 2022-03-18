function ARCMA57_Sham_λ(nlpstop::NLPStopping; λfact::Float64 = 1.0, kwargs...)
    return TRARC(
        nlpstop;
        c = Combi(
            HessSparse,
            PDataMA57,
            (H, g, x, y, z) -> solve_modelARCDiag_HO(H, g, 
                x,
                y,
                z,
                ho_correction = :Shamanskii_MA57,
                λfact = λfact,
            ),
        ),
        kwargs...,
    )
end
