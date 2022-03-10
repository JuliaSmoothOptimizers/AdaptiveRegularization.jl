export TRLDLt_HO_vs_Nwt_λ

function TRLDLt_HO_vs_Nwt_λ(nlp 	   :: AbstractNLPModel,
              	       		nlpstop :: NLPStopping;
				   	        corr_ho :: Symbol = :Shamanskii,
							nwt_res_fact = 0.25,
							λfact = 1.0, 
				   	        kwargs...
               	   	        )

	T = eltype(nlp.meta.x0)

	return TRARC(nlp,
				  nlpstop;
				  TR = TrustRegion(T(10.0)),
				  c = Combi(hessian_dense, PDataLDLt{eltype(nlp.meta.x0)}, (x, y, z) -> solve_modelTRDiag_HO_vs_Nwt_λ(x, y , z; ho_correction = corr_ho, nwt_res_fact = nwt_res_fact, λfact = λfact), preprocessLDLt, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
