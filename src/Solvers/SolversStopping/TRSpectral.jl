export TRSpectral

function TRSpectral(nlp 		:: AbstractNLPModel,
               		nlpstop 	:: NLPStopping;
					kwargs...
               		)

	return TRARC2(nlp,
				  nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi(hessian_dense,PDataSpectral,solve_modelTRDiag,preprocessSpectral,decreaseFact,Tparam()),
				  kwargs...
				  )
end
