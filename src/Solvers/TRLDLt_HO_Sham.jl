export TRLDLt_HO_Sham

function TRLDLt_HO_Sham(nlp 	   :: AbstractNLPModel,
              	   		nlpstop :: NLPStopping;
				   		kwargs...
               	   		)

	return TRARC(nlp,
				  nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi{eltype(nlp.meta.x0)}(hessian_dense, PDataLDLt{eltype(nlp.meta.x0)}, (x, y, z) -> solve_modelTRDiag_HO(x, y, z, ho_correction = :Shamanskii), preprocessLDLt, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end
