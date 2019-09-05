export ST_ARCsparse

function ST_ARCsparse(nlp 		:: AbstractNLPModel,
              	nlpstop 	:: NLPStopping;
				kwargs...
               		)

	return TRARC2(nlp,
				  nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi(hessian_sparse, PDataST{eltype(nlp.meta.x0)}, solve_modelST_ARC, preprocessST, decreaseGen, TparamsST{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
