export TRKsparse

function TRKsparse(nlp 		:: AbstractNLPModel,
               nlpstop 	:: NLPStopping;
							 kwargs...
               )

	return TRARC2(nlp,
				  nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi(hessian_sparse,PDataK,solve_modelKTR,preprocessKTR,decreaseKTR,TparamsKTR()),
				  kwargs...
				  )
end
