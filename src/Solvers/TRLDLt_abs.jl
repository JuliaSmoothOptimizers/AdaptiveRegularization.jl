function TRLDLt_abs(nlpstop 	:: NLPStopping;
					kwargs...
               			)
						   T = eltype(nlpstop.pb.meta.x0)
	return TRARC(nlpstop;
				  TR = TrustRegion(10.0),
				  c = Combi(hessian_dense, PDataLDLt{T}, solve_modelTRDiagAbs, preprocessLDLt, decreaseFact, Tparam{T}()),
				  kwargs...
				  )
end
