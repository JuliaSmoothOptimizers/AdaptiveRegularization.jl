export ST_ARCOp

function ST_ARCOp(nlp 		:: AbstractNLPModel,
              	nlpstop 	:: NLPStopping;
				kwargs...
               		)

	return TRARC2(nlp,
				  nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi(hessian_operator,PDataST,solve_modelST_ARC,preprocessST,decreaseGen,TparamsST()),
				  kwargs...
				  )
end
