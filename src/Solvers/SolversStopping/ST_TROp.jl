export ST_TROp

function ST_TROp(nlp 		:: AbstractNLPModel,
              	nlpstop 	:: NLPStopping;
				kwargs...
               		)

	return TRARC2(nlp,
				  nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi(hessian_operator,PDataST,solve_modelST_TR,preprocessST,decreaseGen,TparamsST()),
				  kwargs...
				  )
end
