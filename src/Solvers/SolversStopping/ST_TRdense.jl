export ST_TRdense

function ST_TRdense(nlp 		:: AbstractNLPModel,
              	nlpstop 	:: NLPStopping;
				kwargs...
               		)

	return TRARC2(nlp,
				  nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi(hessian_dense,PDataST,solve_modelST_TR,preprocessST,decreaseGen,TparamsST()),
				  kwargs...
				  )
end
