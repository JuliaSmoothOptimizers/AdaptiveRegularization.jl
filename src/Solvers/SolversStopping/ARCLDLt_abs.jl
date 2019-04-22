export ARCLDLt_abs

function ARCLDLt_abs(nlp 		:: AbstractNLPModel,
               		nlpstop 	:: NLPStopping;
					kwargs...
               		)

	return TRARC2(nlp,
				  nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi(hessian_dense,PDataLDLt,solve_modelARCDiagAbs,preprocessLDLt,decreaseFact,Tparam()),
				  kwargs...
				  )
end
