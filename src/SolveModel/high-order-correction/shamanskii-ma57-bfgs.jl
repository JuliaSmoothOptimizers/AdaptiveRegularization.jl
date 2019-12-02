export Shamanskii_MA57_BFGS

function Shamanskii_MA57_BFGS(nlp_stop,
					     	  PData :: PDataFact,
					     	  dₙ :: Vector,
						 	  gt :: Vector)
	# printstyled("dans Shamanskii_MA57_BFGS ↓ \n", color = :blue)
	nlp_at_x = nlp_stop.current_state
	∇fₖ = nlp_at_x.gx
	x = nlp_at_x.x
	xdemi = x .+ dₙ
	∇fₖ₊₁ = grad(nlp_stop.pb, xdemi)
	gtemp = copy(∇fₖ₊₁)
	Bₖ = nlp_at_x.Hx
	# @show Bₖ
	# @show PData.L
	# @show PData.D
	# @show PData.s
	# @show PData.pp
	# L = convert(SparseMatrixCSC{Float64, Int64}, PData.L)
	# D = convert(SparseMatrixCSC{Float64, Int64}, PData.D)
	# S = sparse(Diagonal(PData.s))
	# P = SparseMatrixCSC(I, length(gtemp), length(gtemp))[PData.pp, :]
	# Bₖ  = Matrix(inv(Matrix(S)) * P * L * D * L' * P' * inv(Matrix(S)))
	ϵ2 =  sqrt(eps(eltype(gt)))
	Γ = max.(abs.(PData.Δ), ϵ2)
	# αₖ = 1.0
	sₖ = dₙ

	yₖ = ∇fₖ₊₁ .- ∇fₖ
	Yₖ = yₖ * yₖ'
	yts = yₖ' * sₖ
	Yₖ = Yₖ ./ yts

	Bs = Bₖ * sₖ
	BssBₖ = Bs * Bs'
	sBs = Bs' * sₖ
	BssBₖ = (1.0 ./ sBs) .* BssBₖ

	Bₖ₊₁ = Bₖ .+ Yₖ .- BssBₖ
	B57 = convert(SparseMatrixCSC{Float64, Int64}, Bₖ₊₁)
	M = Ma57(B57)
	ma57_factorize(M)
	(L2, D2, s2, pp2) = ma57_get_factors(M)
	# @show L2
	# @show D2
	# @show s2
	# @show pp2
	L2 = SparseMatrixCSC{Float64,Int64}(L2)
	vD1 = diag(D2)       # create internal representation for block diagonal D
    vD2 = diag(D2, 1)     #
    vQ1 = ones(length(vD1))       # vector representation of orthogonal matrix Q
    vQ2 = zeros(length(vD2))      #
    vQ2m = zeros(length(vD2))     #
    veig = copy(vD1)      # vector of eigenvalues of D, initialized to diagonal of D
                          # if D diagonal, nothing more will be computed

    i = 1;
    while i < length(vD1)
        if vD2[i] == 0.0
            i += 1
        else
            mA = [vD1[i] vD2[i]; vD2[i] vD1[i + 1]] #  2X2 submatrix
            # DiagmA, Qma = eig(mA)                   #  spectral decomposition of mA
            X = eigen(mA)
            DiagmA = X.values
            Qma = X.vectors
            veig[i] = DiagmA[1]
            vQ1[i] = Qma[1, 1]
            vQ2[i] = Qma[1, 2]
            vQ2m[i] = Qma[2, 1]
            vQ1[i + 1] = Qma[2, 2]
            veig[i + 1] = DiagmA[2]
            i += 2
        end
    end

    Q2 = sparse(SparseArrays.spdiagm(0 => vQ1, 1 => vQ2m, -1 => vQ2))
	dt = L2 \ gtemp[pp2]
    dc = L2' \ (Q2 * (Q2' * dt ./ Γ))
	d = -dc[invperm(pp2)]

	dₕₒ = dₙ .+ d

	return dₕₒ
end
