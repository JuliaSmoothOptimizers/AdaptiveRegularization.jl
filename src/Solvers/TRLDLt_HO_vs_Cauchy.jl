export TRLDLt_HO_vs_Cauchy

function TRLDLt_HO_vs_Cauchy(nlp 	   :: AbstractNLPModel,
              	   nlpstop :: NLPStopping;
				   corr_ho :: Symbol = :Shamanskii,
				   kwargs...
               	   )

	T = eltype(nlp.meta.x0)
	# printstyled("yo! on est là \n", color = :red)

	return TRARC(nlp,
				  nlpstop;
				  TR = TrustRegion(T(10.0)),
				  c = Combi(hessian_dense, PDataLDLt{eltype(nlp.meta.x0)},  (x, y, z) -> solve_modelTRDiag_HO_vs_Cauchy(x, y , z, ho_correction = corr_ho), preprocessLDLt, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
