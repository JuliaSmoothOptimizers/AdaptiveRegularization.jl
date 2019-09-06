export Shamanskii

function Shamanskii(nlp_stop,
					PData :: PDataFact,
					dₙ :: Vector,
					gt :: Vector)
	# printstyled("on est dans Shamanskii \n", color = :bold)
	nlp_at_x = nlp_stop.current_state
	x = nlp_at_x.x
	# @show x
	# @show dₙ
	xdemi = x + dₙ
	# @show xdemi
	gtemp = grad(nlp_stop.pb, xdemi)
	# @show gtemp

	ϵ2 =  sqrt(eps(eltype(gt)))
	# @show ϵ2
   	Γ = max.(abs.(PData.Δ), ϵ2)
	# @show Γ

	d̃ = PData.L \ gtemp[PData.pp]
	# @show d̃
    d̂ = PData.L' \ (PData.Q * (PData.Q' * d̃ ./ Γ))
	# @show d̂
	d = -d̂[invperm(PData.pp)]
	# @show d
	# printstyled("on sort de Shamanskii \n", color = :green)
	return d, xdemi
end
