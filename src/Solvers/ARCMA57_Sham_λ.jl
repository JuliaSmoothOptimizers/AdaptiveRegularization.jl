function ARCMA57_Sham_λ(nlpstop::NLPStopping; λfact::Float64 = 1.0, kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        c = Combi(
            HessSparse,
            PDataMA57,
            (x, y, z) -> solve_modelARCDiag_HO(
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
