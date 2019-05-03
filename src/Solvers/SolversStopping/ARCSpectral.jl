export ARCSpectral

function ARCSpectral(nlp 		:: AbstractNLPModel,
               		 nlpstop 	:: NLPStopping;
					 kwargs...
               		 )

	return TRARC2(nlp,
				  nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi(ARCTR.hessian_dense,ARCTR.PDataSpectral{eltype(nlp.meta.x0)},ARCTR.solve_modelARCDiag,ARCTR.preprocessSpectral,ARCTR.decreaseFact,ARCTR.Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
