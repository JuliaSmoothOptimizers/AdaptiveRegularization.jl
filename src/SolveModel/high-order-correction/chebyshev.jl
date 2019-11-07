export Chebyshev

function Chebyshev(nlp_stop,
				   PData :: PDataFact,
				   dₙ :: Vector,
				   gt :: Vector)
	# printstyled("on est dans Chebyshev  \n", color = :yellow)
	T = eltype(gt)
	# @show T
	nlp_at_x = nlp_stop.current_state
	x = nlp_at_x.x
	# xdemi = NaN * rand(length(dₙ))
	# @show x
	# @show ∇f³xuv
	# @show typeof(nlp_stop.pb)
	Tdndn = ∇f³xuv(nlp_stop.pb, x, dₙ, dₙ)
	# printstyled("on a Tdndn  \n", color = :yellow)

	ϵ2 =  sqrt(eps(eltype(gt)))
   	Γ = max.(abs.(PData.Δ), ϵ2)

	d̃ = PData.L \ Tdndn[PData.pp]
    d̂ = PData.L' \ (PData.Q * (PData.Q' * d̃ ./ Γ))
    d = -d̂[invperm(PData.pp)]

	dC = dₙ + T(0.5) * d
	return dC
end
