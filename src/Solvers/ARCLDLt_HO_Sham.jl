function ARCLDLt_HO_Sham(nlpstop 	:: NLPStopping;
						λfact = 100.0,
						kwargs...
               			)

						   T = eltype(nlpstop.pb.meta.x0)

	return TRARC(nlpstop;
				  TR = TrustRegion(T(10.0)),
				  c = Combi(hessian_dense, PDataLDLt{T}, (x, y, z) -> solve_modelARCDiag_HO(x, y, z, ho_correction = :Shamanskii, λfact = λfact), preprocessLDLt, decreaseFact, Tparam{T}()),
				  kwargs...
				  )
end
