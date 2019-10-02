export TRMA57_Sham_2

function TRMA57_Sham_2(nlp 		:: AbstractNLPModel,
              		   nlpstop 	:: NLPStopping;
					   kwargs...
               		   )

	return TRARC(nlp, nlpstop; TR = TrustRegion(10.0),
				  c = Combi(hessian_sparse, PDataMA57{eltype(nlp.meta.x0)}, (x, y, z) -> solve_modelTRDiag_HO(x, y, z, ho_correction = :Shamanskii_MA57), preprocessMA57_2, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
