export TRMA57_Sham

function TRMA57_Sham(nlp 		:: AbstractNLPModel,
              		 nlpstop 	:: NLPStopping;
					 kwargs...
               		 )

	# printstyled("on est dans TRMA57_Sham \n", color = :blue)

	return TRARC2(nlp, nlpstop; TR = TrustRegion(10.0),
				  c = Combi(hessian_sparse, PDataMA57{eltype(nlp.meta.x0)}, solve_modelTRDiag_HO, preprocessMA57, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
				  kwargs...
				  )
end

# function TRMA57_Sham(nlp 		:: AbstractNLPModel,
#               		 nlpstop 	:: NLPStopping;
# 					 kwargs...
#                		 )
#
# 	return TRARC2(nlp, nlpstop; TR = TrustRegion(10.0),
# 				  c = Combi(hessian_sparse, PDataMA57{eltype(nlp.meta.x0)}, solve_modelTRDiag, preprocessMA57, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
# 				  kwargs...
# 				  )
# end
