export TRSpectral

function TRSpectral(nlp 		:: AbstractNLPModel,
               		nlpstop 	:: NLPStopping;
					kwargs...
               		)

	return TRARC(nlp,
				  nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi(hessian_dense, PDataSpectral{eltype(nlp.meta.x0)}, solve_modelTRDiag, preprocessSpectral, decreaseFact,Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
