export TRLDLt_HO

function TRLDLt_HO(nlp 	   :: AbstractNLPModel,
              	   nlpstop :: NLPStopping;
				   corr_ho :: Symbol = :Shamanskii,
				   kwargs...
               	   )

	return TRARC(nlp,
				  nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi(hessian_dense,PDataLDLt{eltype(nlp.meta.x0)},  (x, y, z) -> solve_modelTRDiag_HO(x, y , z, ho_correction = corr_ho, fact = 2.0), preprocessLDLt, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
