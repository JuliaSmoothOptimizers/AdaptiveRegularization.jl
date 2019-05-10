export Shamanskii

function Shamanskii(nlp_stop,
					PData :: PDataFact,
					dₙ :: Vector,
					gt :: Vector)
	nlp_at_x = nlp_stop.current_state
	x = nlp_at_x.x
	xdemi = x + dₙ
	gtemp = grad(nlp_stop.pb, xdemi)

	ϵ2 =  sqrt(eps(eltype(gt)))
   	Γ = max.(abs.(PData.Δ), ϵ2)

	d̃ = PData.L \ gtemp[PData.pp]
    d̂ = PData.L' \ (PData.Q * (PData.Q' * d̃ ./ Γ))
    d = -d̂[invperm(PData.pp)]
	return d, xdemi
end
