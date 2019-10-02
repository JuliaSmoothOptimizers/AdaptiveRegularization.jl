export ARCLDLt_abs

function ARCLDLt_abs(nlp 		:: AbstractNLPModel,
               		nlpstop 	:: NLPStopping;
					kwargs...
               		)

	return TRARC(nlp,
				  nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi(hessian_dense, PDataLDLt{eltype(nlp.meta.x0)}, solve_modelARCDiagAbs, preprocessLDLt, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
