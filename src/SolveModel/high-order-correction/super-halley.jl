export SuperHalley

function SuperHalley(nlp_stop, PData::PDataFact, dₙ::Vector, gt::Vector)

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

    # hess_plus_tenseur = (nlp_at_x.Hx + ∇³fdₙ)
    if PData.λ .== 0.0
        hess_plus_tenseur = (nlp_at_x.Hx .+ ∇³fdₙ)
    else
        hess_plus_tenseur = (nlp_at_x.Hx .+ ∇³fdₙ) .+ (PData.λ .* Matrix(1.0I, n, n))
    end

    if isposdef(hess_plus_tenseur)
        (L, D, pp, ρ, ncomp) = ldlt_symm(hess_plus_tenseur, 'r')
    else
        # println("on est ici")
        return dₙ
    end

    # try
    #    (L, D, pp, ρ, ncomp) = ldlt_symm(hess_plus_tenseur, 'r')
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
    Q = X.vectors

    ϵ = sqrt(eps(eltype(gt)))
    Γ = max.(abs.(Δ), ϵ)

    d̃ = L \ nlp_at_x.gx[pp]
    d̂ = L' \ (Q * (Q' * d̃ ./ Γ))
    dV = -d̂[invperm(pp)]

    grad_plus_tuv = -(nlp_at_x.gx + 0.5 .* ∇f³xuv(nlp_stop.pb, x, dₙ, dV))
    grad_plus_tuv = -grad_plus_tuv

    ϵ2 = sqrt(eps(eltype(gt)))
    Γ2 = max.(abs.(PData.Δ), ϵ2)

    d̃ = PData.L \ grad_plus_tuv[pp]
    d̂ = PData.L' \ (PData.Q * (PData.Q' * d̃ ./ Γ2))
    dSH = -d̂[invperm(PData.pp)]

    return dSH
end
