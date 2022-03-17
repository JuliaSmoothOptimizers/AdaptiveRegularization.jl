function TRMA97(nlpstop::NLPStopping; kwargs...)
    T = eltype(nlpstop.pb.meta.x0)
    return TRARC(
        nlpstop;
        c = Combi(
            HessSparse,
            PDataMA97,
            solve_modelTRDiag,
        ),
        kwargs...,
    )
end
