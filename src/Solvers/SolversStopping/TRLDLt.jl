export TRLDLt

function TRLDLt(nlp 		:: AbstractNLPModel,
              	nlpstop 	:: NLPStopping;
				kwargs...
               		)

	return TRARC2(nlp,
				  nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi(hessian_dense,PDataLDLt,solve_modelTRDiag,preprocessLDLt,decreaseFact,Tparam()),
				  kwargs...
				  )
end
