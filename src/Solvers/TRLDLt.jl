export TRLDLt

function TRLDLt(nlp 		:: AbstractNLPModel,
              	nlpstop 	:: NLPStopping;
				kwargs...
               		)

	T = eltype(nlp.meta.x0)

	return TRARC(nlp,
				  nlpstop;
				  TR = TrustRegion(T(10.0)),
				  c = Combi(hessian_dense, PDataLDLt{eltype(nlp.meta.x0)}, solve_modelTRDiag, preprocessLDLt, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
