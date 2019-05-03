export eigen2x2
function eigen2x2(A)
	Tp = eltype(A)
	a = A[1, 1]; b = A[1, 2]; c = A[2, 1]; d = A[2, 2];
	T = tr(A); D = det(A);

	L₁ = T / 2 + sqrt((T^2 / 4) - D)
	L₂ = T / 2 - sqrt((T^2 / 4) - D)
	values = [L₁, L₂]

	global vectors = []

	if !iszero(c)
		push!(vectors, [L₁ - d, c])
		push!(vectors, [L₂ - d, c])
	elseif !iszero(b)
		push!(vectors, [b, L₁ - a])
		push!(vectors, [b, L₂ - a])
	elseif iszero(b) && iszero(c)
		push!(vectors, [one(Tp), zero(Tp)])
		push!(vectors, [zero(Tp), one(Tp)])
	end

	return values, vectors
end
