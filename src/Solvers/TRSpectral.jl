function TRSpectral(nlpstop 	:: NLPStopping;
					kwargs...
               		)
					   T = eltype(nlpstop.pb.meta.x0)
	return TRARC(nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi(hessian_dense, PDataSpectral{T}, solve_modelTRDiag, preprocessSpectral, decreaseFact,Tparam{T}()),
				  kwargs...
				  )
end
