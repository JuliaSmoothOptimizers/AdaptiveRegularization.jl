export TRMA57_2

function TRMA57_2(nlp 		:: AbstractNLPModel,
              	  nlpstop 	:: NLPStopping;
				  kwargs...
               	  )

	return TRARC(nlp, nlpstop; TR = TrustRegion(10.0),
				  c = Combi(hessian_sparse_no_tril, PDataMA57{eltype(nlp.meta.x0)}, solve_modelTRDiag, preprocessMA57, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
