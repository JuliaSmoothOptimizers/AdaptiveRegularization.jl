function ST_ARCdense(nlpstop 	:: NLPStopping;
				kwargs...
               		)
					   T = eltype(nlpstop.pb.meta.x0)
	return TRARC(nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi(hessian_dense, PDataST{T}, solve_modelST_ARC, preprocessST, decreaseGen, TparamsST{T}()),
				  kwargs...
				  )
end
