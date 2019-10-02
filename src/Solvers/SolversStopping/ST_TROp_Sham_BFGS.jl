export ST_TROp_Sham_BFGS

function ST_TROp_Sham_BFGS(nlp 		:: AbstractNLPModel,
              	nlpstop 	:: NLPStopping;
				kwargs...
               		)

	return TRARC2(nlp,
				  nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi(hessian_operator, PDataST{eltype(nlp.meta.x0)}, solve_modelST_TR_Sham_BFGS, preprocessST, decreaseGen, TparamsST{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
