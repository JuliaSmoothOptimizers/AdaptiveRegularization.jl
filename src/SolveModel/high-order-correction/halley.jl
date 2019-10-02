export Halley

function Halley(nlp_stop,
				   PData :: PDataFact,
				   dₙ :: Vector,
				   gt :: Vector)

	nlp_at_x = nlp_stop.current_state
	x = nlp_at_x.x

	n = length(x)

	e₁ = zeros(n)
    e₁[1] = 1.0
    ∇³fdₙ = ∇f³xuv(nlp_stop.pb, x, dₙ, e₁)
    for i = 2:n
        eᵢ = zeros(n)
        eᵢ[i] = 1.0
        ∇³fdₙ = hcat(∇³fdₙ, ∇f³xuv(nlp_stop.pb, x, dₙ, eᵢ))
    end

	hess_plus_tenseur = (nlp_at_x.Hx + 0.5   .* ∇³fdₙ)

	try
       (L, D, pp, ρ, ncomp) = ldlt_symm(hess_plus_tenseur, 'r')
	catch
 		println("*******   Problem in LDLt Halley")
       	#res = PDataLDLt()
       	#res.OK = false
       	#return res
       	d = NaN * zeros(n)
       	return d, xdemi
   	end

	if true in isnan.(D)
   	   println("*******   Problem in D from LDLt: NaN")
	   println(" cond (H) = $(cond(H))")
	   #res = PDataLDLt()
	   #res.OK = false
	   #return res
	   d = NaN * zeros(n)
	   return d, xdemi
   	end

	D = Hermitian(D)
   	X = eigen(D)
  	Δ = X.values
  	Q =  X.vectors

	ϵ2 =  sqrt(eps(eltype(gt)))
   	Γ = max.(abs.(Δ), ϵ2)

	d̃ = L \ nlp_at_x.gx[pp]
    d̂ = L' \ (Q * (Q' * d̃ ./ Γ))
    dH = -d̂[invperm(pp)]

	return dH
end
