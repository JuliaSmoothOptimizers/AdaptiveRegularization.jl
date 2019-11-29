export Halley

function Halley(nlp_stop,
				   PData :: PDataFact,
				   dₙ :: Vector,
				   gt :: Vector)
	# println("on est dans Halley")
	nlp_at_x = nlp_stop.current_state
	x = nlp_at_x.x
	# println("on a nlp_at_x et x")

	n = length(x)
	# @show n

	e₁ = zeros(n)
	# @show e₁
    e₁[1] = 1.0
	# @show e₁
    ∇³fdₙ = ARCTR.∇f³xuv(nlp_stop.pb, x, dₙ, e₁)
	# @show ∇³fdₙ
	# println("avant de rentrer dans le for")
	# @show n
    for i = 2:n
		# @show i
        eᵢ = zeros(n)
        eᵢ[i] = 1.0
        ∇³fdₙ = hcat(∇³fdₙ, ∇f³xuv(nlp_stop.pb, x, dₙ, eᵢ))
    end

	# @show PData.λ
	# @show PData.λ .== 0.0
	if PData.λ .== 0.0
		hess_plus_tenseur = (nlp_at_x.Hx .+ 0.5 .* ∇³fdₙ)
	else
		hess_plus_tenseur = (nlp_at_x.Hx .+ 0.5 .* ∇³fdₙ) .+ (PData.λ .* Matrix(1.0I, n, n))
	end
	# hess_plus_tenseur = (nlp_at_x.Hx .+ 0.5 .* ∇³fdₙ)
	# @show isposdef(hess_plus_tenseur)
	if isposdef(hess_plus_tenseur)
        (L, D, pp, ρ, ncomp) = ldlt_symm(hess_plus_tenseur, 'r')
   	else
		# println("on est ici")
	   	return dₙ
   	end
	# catch
 	# 	println("*******   Problem in LDLt Halley")
    #    	#res = PDataLDLt()
    #    	#res.OK = false
    #    	#return res
    #    	d = NaN * zeros(n)
    #    	return d, xdemi
   	# end

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
