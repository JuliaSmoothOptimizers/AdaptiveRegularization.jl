export TRKsparse

function TRKsparse(nlp 		:: AbstractNLPModel,
               nlpstop 	:: NLPStopping;
							 kwargs...
               )

	return TRARC(nlp,
				  nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi(hessian_sparse, PDataK{eltype(nlp.meta.x0)}, solve_modelKTR, preprocessKTR, decreaseKTR, TparamsKTR{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
