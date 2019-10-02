export ARCqKOp1


function ARCqKOp1(nlp 		:: AbstractNLPModel,
              	nlpstop 	:: NLPStopping;
				kwargs...
               		)

    shifts = 10.0.^(collect(-15.0:1.0:15.0))

    return TRARC2(nlp,
		  nlpstop;
		  TR = TrustRegion(10.0),
		  c = Combi(hessian_operator, PDataK{eltype(nlp.meta.x0)}, solve_modelKARC, preprocessKARC, decreaseKARC, TparamsKARC{eltype(nlp.meta.x0)}(shifts)),
		  kwargs...
		  )
end
