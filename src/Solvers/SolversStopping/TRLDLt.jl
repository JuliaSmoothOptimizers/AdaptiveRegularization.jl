export TRLDLt

function TRLDLt(nlp 		:: AbstractNLPModel,
              	nlpstop 	:: NLPStopping;
				kwargs...
               		)

	return TRARC2(nlp,
				  nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi(hessian_dense,PDataLDLt{eltype(nlp.meta.x0)},solve_modelTRDiag,preprocessLDLt,decreaseFact,Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
