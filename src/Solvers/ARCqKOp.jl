export ARCqKOp

function ARCqKOp(nlp 		:: AbstractNLPModel,
              	nlpstop 	:: NLPStopping;
				kwargs...
               		)

	return TRARC(nlp,
				  nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi(hessian_operator, PDataK{eltype(nlp.meta.x0)}, solve_modelKARC, preprocessKARC, decreaseKARC, TparamsKARC{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
