export TRMA57_Sham

function TRMA57_Sham(nlp 		:: AbstractNLPModel,
              		 nlpstop 	:: NLPStopping;
					 kwargs...
               		 )

	return TRARC2(nlp, nlpstop; TR = TrustRegion(10.0),
				  c = Combi(hessian_sparse, PDataMA57{eltype(nlp.meta.x0)}, (x, y, z) -> solve_modelTRDiag_HO(x, y, z, ho_correction = :Shamanskii_MA57), preprocessMA57, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
