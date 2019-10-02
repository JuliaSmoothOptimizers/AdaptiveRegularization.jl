export ARCMA57

function ARCMA57(nlp 		:: AbstractNLPModel,
              	 nlpstop 	:: NLPStopping;
				 kwargs...
               	 )

	return TRARC(nlp, nlpstop; TR = TrustRegion(10.0),
				  c = Combi(hessian_sparse, PDataMA57{eltype(nlp.meta.x0)}, solve_modelARCDiag, preprocessMA57, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
