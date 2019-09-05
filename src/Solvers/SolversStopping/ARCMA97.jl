export ARCMA97

function ARCMA97(nlp 		:: AbstractNLPModel,
              	 nlpstop 	:: NLPStopping;
				 kwargs...
               	 )

	return TRARC2(nlp, nlpstop; TR = TrustRegion(10.0),
				  c = Combi(hessian_dense, PDataMA97{eltype(nlp.meta.x0)}, solve_modelARCDiag, preprocessMA97, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
