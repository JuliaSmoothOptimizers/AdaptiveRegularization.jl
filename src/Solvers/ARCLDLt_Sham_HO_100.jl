export ARCLDLt_HO_Sham_100

function ARCLDLt_HO_Sham_100(nlp 		:: AbstractNLPModel,
              				 nlpstop 	:: NLPStopping;
							 kwargs...
               				 )

	T = eltype(nlp.meta.x0)

	return TRARC(nlp,
				  nlpstop;
				  TR = TrustRegion(T(10.0)),
				  c = Combi(hessian_dense, PDataLDLt{eltype(nlp.meta.x0)}, (x, y, z) -> solve_modelARCDiag_HO(x, y, z, λfact = 100.0), preprocessLDLt, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
