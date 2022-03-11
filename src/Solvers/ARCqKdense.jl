function ARCqKdense(nlpstop 	:: NLPStopping;
				kwargs...
               		)
					   T = eltype(nlpstop.pb.meta.x0)
	shifts = 10.0.^(collect(-15.0:1.0:15.0))

	return TRARC(nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi(hessian_dense, PDataK{T}, solve_modelKARC, preprocessKARC, decreaseKARC, TparamsKARC{T}(shifts)),
				  kwargs...
				  )
end
