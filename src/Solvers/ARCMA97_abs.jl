export ARCMA97_abs

function ARCMA97_abs(nlp 		:: AbstractNLPModel,
              	 	 nlpstop 	:: NLPStopping;
				 	 kwargs...
               	 	 )

	return TRARC(nlp, nlpstop; TR = TrustRegion(10.0),
				  c = Combi(hessian_dense, PDataMA97{eltype(nlp.meta.x0)}, solve_modelARCDiagAbs, preprocessMA97, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
