export TRMA57_abs

function TRMA57_abs(nlp 		:: AbstractNLPModel,
              		nlpstop 	:: NLPStopping;
					kwargs...
               		)

	return TRARC2(nlp, nlpstop; TR = TrustRegion(10.0),
				  c = Combi(hessian_sparse, PDataMA57{eltype(nlp.meta.x0)}, solve_modelTRDiagAbs, preprocessMA57, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
