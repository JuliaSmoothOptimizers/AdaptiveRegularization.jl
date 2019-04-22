export TRKdense

function TRKdense(nlp 		:: AbstractNLPModel,
               nlpstop 	:: NLPStopping;
							 kwargs...
               )

	return TRARC2(nlp,
				  nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi(hessian_dense,PDataK,solve_modelKTR,preprocessKTR,decreaseKTR,TparamsKTR()),
				  kwargs...
				  )
end
