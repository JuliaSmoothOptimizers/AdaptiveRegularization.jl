function TRKdense(nlpstop 	:: NLPStopping;
							 kwargs...
               )
			   T = eltype(nlpstop.pb.meta.x0)
	return TRARC(nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi(hessian_dense, PDataK{T}, solve_modelKTR, preprocessKTR, decreaseKTR, TparamsKTR{T}()),
				  kwargs...
				  )
end
