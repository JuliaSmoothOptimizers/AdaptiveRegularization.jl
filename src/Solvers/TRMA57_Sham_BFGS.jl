export TRMA57_Sham_BFGS

function TRMA57_Sham_BFGS(nlp 		:: AbstractNLPModel,
              		 	  nlpstop 	:: NLPStopping;
					 	  kwargs...
               		 	  )

	return TRARC(nlp, nlpstop; TR = TrustRegion(10.0),
				  c = Combi(hessian_sparse, PDataMA57{eltype(nlp.meta.x0)}, (x, y, z) -> solve_modelTRDiag_HO(x, y, z, ho_correction = :Shamanskii_MA57_BFGS), preprocessMA57, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
