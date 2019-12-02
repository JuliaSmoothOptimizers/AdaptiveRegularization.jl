export ARCLDLt

function ARCLDLt(nlp 		:: AbstractNLPModel,
              	nlpstop 	:: NLPStopping;
				kwargs...
               		)

	T = eltype(nlp.meta.x0)

	return TRARC(nlp,
				  nlpstop;
				  TR = TrustRegion(T(10.0)),
				  c = Combi(hessian_dense, PDataLDLt{eltype(nlp.meta.x0)}, solve_modelARCDiag, preprocessLDLt, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
