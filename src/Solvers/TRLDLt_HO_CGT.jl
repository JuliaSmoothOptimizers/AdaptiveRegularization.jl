export TRLDLt_HO_CGT

function TRLDLt_HO_CGT(nlp 	   :: AbstractNLPModel,
              	       nlpstop :: NLPStopping;
				   	   corr_ho :: Symbol = :Shamanskii,
					   κmdc = 1e-08,
				   	   kwargs...
               	   	   )

	T = eltype(nlp.meta.x0)
	κmdc = T(κmdc)

	return TRARC(nlp,
				  nlpstop;
				  TR = TrustRegion(T(10.0)),
				  c = Combi(hessian_dense, PDataLDLt{eltype(nlp.meta.x0)},  (x, y, z) -> solve_modelTRDiag_HO_CGT(x, y , z, ho_correction = corr_ho, κmdc = κmdc), preprocessLDLt, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
