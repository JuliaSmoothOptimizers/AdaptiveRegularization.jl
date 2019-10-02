export Shamanskii_BFGS

function Shamanskii_BFGS(nlp_stop,
					     PData  :: PDataFact,
					     dₙ 	:: Vector,
					     gt 	:: Vector)
	nlp_at_x = nlp_stop.current_state
	∇fₖ = nlp_at_x.gx
	x = nlp_at_x.x
	Bₖ = Symmetric(hess(nlp_stop.pb, x), :L)

	sₖ = dₙ
	xtemp = x .+ sₖ
	∇fₖ₊₁ = grad(nlp_stop.pb, xtemp)
	yₖ = ∇fₖ₊₁ .- ∇fₖ
	Yₖ = yₖ * yₖ'
	yts = yₖ' * sₖ
	Yₖ = Yₖ ./ yts

	Bs = Bₖ * sₖ
	BssBₖ = Bs * Bs'
	sBs = Bs' * sₖ
	BssBₖ = (1.0 ./ sBs) .* BssBₖ

	Bₖ₊₁ = Bₖ .+ Yₖ .- BssBₖ
	LDLTB = ldl(Bₖ₊₁)
	dₛ = -(LDLTB \ grad(nlp_stop.pb, xtemp))
	xt = xtemp + dₛ

	dₕₒ = dₙ .+ dₛ

	return dₕₒ
end
