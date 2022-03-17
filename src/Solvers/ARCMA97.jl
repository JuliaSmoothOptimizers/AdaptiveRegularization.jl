function ARCMA97(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        c = Combi(
            HessDense,
            PDataMA97,
            solve_modelARCDiag,
        ),
        kwargs...,
    )
end
