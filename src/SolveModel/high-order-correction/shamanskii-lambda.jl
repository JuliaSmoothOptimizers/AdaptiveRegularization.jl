export Shamanskii_λ

function Shamanskii_λ(nlp_stop,
					  PData :: PDataFact,
					  dₙ :: Vector,
					  gt :: Vector;
					  λnull :: Bool = false)
	nlp_at_x = nlp_stop.current_state
	x = nlp_at_x.x
	n = length(x)
	xdemi = x + dₙ
	gtemp = grad(nlp_stop.pb, xdemi)

	ϵ2 =  sqrt(eps(eltype(gt)))
   	Γ = max.(abs.(PData.Δ), ϵ2)

	if λnull
		d̃ = PData.L \ gtemp[PData.pp]
    	d̂ = PData.L' \ (PData.Q * (PData.Q' * d̃ ./ Γ))
		d = -d̂[invperm(PData.pp)]
	else
		D2 = SymTridiagonal(Hermitian(PData.D + PData.λ * Matrix(1.0I, n, n)))
		Q2 = eigvecs(D2)
		d̃ = PData.L \ gtemp[PData.pp]
    	d̂ = PData.L' \ (Q2 * (Q2' * d̃ ./ Γ))
		d = -d̂[invperm(PData.pp)]
	end

	dₕₒ = dₙ .+ d

	return dₕₒ
end
