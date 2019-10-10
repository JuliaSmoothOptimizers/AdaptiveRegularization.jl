export TRLDLt_HO_Sham_λ

function TRLDLt_HO_Sham_λ(nlp 	   :: AbstractNLPModel,
              	   		nlpstop :: NLPStopping;
				   		kwargs...
               	   		)

	return TRARC(nlp,
				  nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi{eltype(nlp.meta.x0)}(hessian_dense, PDataLDLt{eltype(nlp.meta.x0)}, (x, y, z) -> solve_modelTRDiag_HO_λ(x, y, z, ho_correction = :Shamanskii_λ, λfact = 100.0), preprocessLDLt, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
