export TRKOp

function TRKOp(nlp 		:: AbstractNLPModel,
               nlpstop 	:: NLPStopping;
							 kwargs...
               )

	return TRARC2(nlp,
				  nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi(hessian_operator, PDataK{eltype(nlp.meta.x0)}, solve_modelKTR, preprocessKTR, decreaseKTR, TparamsKTR{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
