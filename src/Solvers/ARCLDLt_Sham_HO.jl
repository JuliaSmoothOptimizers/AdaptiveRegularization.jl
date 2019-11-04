export ARCLDLt_HO_Sham

function ARCLDLt_HO_Sham(nlp 		:: AbstractNLPModel,
              			nlpstop 	:: NLPStopping;
						λfact = 100.0,
						kwargs...
               			)

	T = eltype(nlp.meta.x0)

	return TRARC(nlp,
				  nlpstop;
				  TR = TrustRegion(T(10.0)),
				  c = Combi(hessian_dense, PDataLDLt{eltype(nlp.meta.x0)}, (x, y, z) -> solve_modelARCDiag_HO(x, y, z, ho_correction = :Shamanskii, λfact = λfact), preprocessLDLt, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
