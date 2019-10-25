export ARCLDLt_HO_vs_Nwt

function ARCLDLt_HO_vs_Nwt(nlp 	   :: AbstractNLPModel,
              	   nlpstop :: NLPStopping;
				   corr_ho :: Symbol = :Shamanskii,
				   nwt_res_fact = 0.8,
				   λfact = 100.0,
				   kwargs...
               	   )

	T = eltype(nlp.meta.x0)

	return TRARC(nlp,
				  nlpstop;
				  TR = TrustRegion(T(10.0)),
				  c = Combi(hessian_dense, PDataLDLt{eltype(nlp.meta.x0)},  (x, y, z) -> solve_modelARCDiag_HO_vs_Nwt(x, y, z, λfact = λfact, nwt_res_fact = nwt_res_fact), preprocessLDLt, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
