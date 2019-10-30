export TRMA57_Sham_vs_Nwt

function TRMA57_Sham_vs_Nwt(nlp 		:: AbstractNLPModel,
              		 		nlpstop 	:: NLPStopping;
							corr_ho :: Symbol = :Shamanskii_MA57,
							nwt_res_fact = 0.25,
					 		kwargs...
               		 		)

	return TRARC(nlp, nlpstop; TR = TrustRegion(10.0),
				  c = Combi(hessian_sparse, PDataMA57{eltype(nlp.meta.x0)}, (x, y, z) -> solve_modelTRDiag_HO_vs_Nwt(x, y, z; ho_correction = corr_ho, nwt_res_fact = nwt_res_fact), preprocessMA57, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
