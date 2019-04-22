export ARCqKOp

function ARCqKOp(nlp 		:: AbstractNLPModel,
              	nlpstop 	:: NLPStopping;
				kwargs...
               		)

	return TRARC2(nlp,
				  nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi(hessian_operator,PDataK,solve_modelKARC,preprocessKARC,decreaseKARC,TparamsKARC()),
				  kwargs...
				  )
end
