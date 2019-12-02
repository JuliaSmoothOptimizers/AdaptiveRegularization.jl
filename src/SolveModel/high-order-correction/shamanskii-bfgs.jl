export Shamanskii_BFGS

function Shamanskii_BFGS(nlp_stop,
					     PData  :: PDataFact,
					     dₙ 	:: Vector,
					     gt 	:: Vector)
	# printstyled("on est ici ↓ \n", color = :bold)
	T = eltype(gt)
	nlp_at_x = nlp_stop.current_state
	∇fₖ = nlp_at_x.gx
	x = nlp_at_x.x
	# Bₖ = Symmetric(hess(nlp_stop.pb, x), :L)
	# @show nlp_at_x.Hx
	Bₖ = nlp_at_x.Hx
	# printstyled("on a Bₖ \n", color = :bold)

	sₖ = dₙ
	# printstyled("on a sₖ \n", color = :bold)
	xtemp = x .+ sₖ
	# printstyled("on a xtemp \n", color = :bold)
	∇fₖ₊₁ = grad(nlp_stop.pb, xtemp)
	# printstyled("on a fk+1 \n", color = :bold)
	yₖ = ∇fₖ₊₁ .- ∇fₖ

	Yₖ = yₖ * yₖ'
	yts = yₖ' * sₖ
	Yₖ = Yₖ ./ yts
	# printstyled("on a Yₖ \n", color = :bold)

	Bs = Bₖ * sₖ
	BssBₖ = Bs * Bs'
	sBs = Bs' * sₖ
	BssBₖ = (1.0 ./ sBs) .* BssBₖ

	Bₖ₊₁ = Bₖ .+ Yₖ .- BssBₖ
	# LDLTB = ldl(Bₖ₊₁)
	L, D, pp, ρ, ncomp = ldlt_symm(Bₖ₊₁, 'r')
    # printstyled("on a L, D,... \n", color = :bold)
	D = Hermitian(D)
    X = eigen(D)
	# printstyled("on a X \n", color = :bold)
    Δ = X.values
	Q =  X.vectors
    l_m, = findmin(Δ)
	# printstyled("on a l_m \n", color = :bold)
    ĝ = L \ (∇fₖ₊₁[pp])
	# printstyled("on a g chapeau \n", color = :bold)
    g̃ = Q' * ĝ
	# printstyled("on a g tilde \n", color = :bold)
    n_g = norm(∇fₖ₊₁)
    λ =  max(-l_m, 0.0)
	# printstyled("on a λ \n", color = :bold)
	ϵ2 = sqrt(eps(T))
	Γ = max.(abs.(Δ), ϵ2)

	d̃ = L \ ∇fₖ₊₁[pp]
   	d̂ = L' \ (Q * (Q' * d̃ ./ Γ))
   	dₛ = -d̂[invperm(pp)]
	# xt = xtemp + dₛ

	dₕₒ = dₙ .+ dₛ

	return dₕₒ
end
