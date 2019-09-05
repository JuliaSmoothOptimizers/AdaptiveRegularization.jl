export ST_ARCdense

function ST_ARCdense(nlp 		:: AbstractNLPModel,
              	nlpstop 	:: NLPStopping;
				kwargs...
               		)

	return TRARC2(nlp,
				  nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi(hessian_dense, PDataST{eltype(nlp.meta.x0)}, solve_modelST_ARC, preprocessST, decreaseGen, TparamsST{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
