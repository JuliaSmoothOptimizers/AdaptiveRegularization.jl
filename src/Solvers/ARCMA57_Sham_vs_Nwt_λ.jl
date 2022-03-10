export ARCMA57_Sham_vs_Nwt_λ

function ARCMA57_Sham_vs_Nwt_λ(nlp 		:: AbstractNLPModel,
              	 			   nlpstop 	:: NLPStopping;
							   λfact 	:: Float64 = 1.0,
							   nwt_res_fact = 0.25,
				 			   kwargs...
               	 			   )

	return TRARC(nlp, nlpstop; TR = TrustRegion(10.0),
				  c = Combi(hessian_sparse, PDataMA57{eltype(nlp.meta.x0)}, (x, y, z) -> solve_modelARCDiag_HO_vs_Nwt(x, y, z, λfact = λfact, nwt_res_fact = nwt_res_fact, ho_correction = :Shamanskii_MA57), preprocessMA57, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
