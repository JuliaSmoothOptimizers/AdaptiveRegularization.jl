export TRMA57_Sham_λ

function TRMA57_Sham_λ(nlp 		:: AbstractNLPModel,
              		   nlpstop 	:: NLPStopping;
					   λfact    :: Float64 = 10.0,
					   kwargs...
               		   )

	return TRARC(nlp, nlpstop; TR = TrustRegion(10.0),
				  c = Combi(hessian_sparse, PDataMA57{eltype(nlp.meta.x0)}, (x, y, z) -> solve_modelTRDiag_HO_λ(x, y, z, ho_correction = :Shamanskii_MA57, λfact = λfact), preprocessMA57, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
