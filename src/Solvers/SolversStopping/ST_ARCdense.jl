export ST_ARCdense

function ST_ARCdense(nlp 		:: AbstractNLPModel,
              	nlpstop 	:: NLPStopping;
				kwargs...
               		)

	return TRARC2(nlp,
				  nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi(hessian_dense,PDataST,solve_modelST_ARC,preprocessST,decreaseGen,TparamsST()),
				  kwargs...
				  )
end
