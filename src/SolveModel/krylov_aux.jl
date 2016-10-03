"""Numerically stable symmetric Givens rotation.
Given `a` and `b`, return `(c, s, ρ)` such that

    [ c  s ] [ a ] = [ ρ ]
    [ s -c ] [ b ] = [ 0 ].
"""
function sym_givens(a :: Float64, b :: Float64)
	#
	# Modeled after the corresponding Matlab function by M. A. Saunders and S.-C. Choi.
	# http://www.stanford.edu/group/SOL/dissertations/sou-cheng-choi-thesis.pdf
	# D. Orban, Montreal, May 2015.

  if b == 0.0
    a == 0.0 && (c = 1.0) || (c = sign(a));  # In Julia, sign(0) = 0.
    s = 0.0;
    ρ = abs(a);

  elseif a == 0.0
    c = 0.0;
    s = sign(b);
    ρ = abs(b);

  elseif abs(b) > abs(a)
    t = a / b;
    s = sign(b) / sqrt(1.0 + t * t);
    c = s * t;
    ρ = b / s;  # Computationally better than d = a / c since |c| <= |s|.

  else
    t = b / a;
    c = sign(a) / sqrt(1.0 + t * t);
    s = c * t;
    ρ = a / c;  # Computationally better than d = b / s since |s| <= |c|
  end

  return (c, s, ρ)
end


"""Find the real roots of the quadratic

    q(x) = q₂ x² + q₁ x + q₀,

where q₂, q₁ and q₀ are real. Care is taken to avoid numerical
cancellation. Optionally, `nitref` steps of iterative refinement
may be performed to improve accuracy. By default, `nitref=1`.
"""
function roots_quadratic(q₂ :: Float64, q₁ :: Float64, q₀ :: Float64;
                         nitref :: Int=1)
  # Case where q(x) is linear.
  if q₂ == 0.0
    if q₁ == 0.0
      q₀ == 0.0 && return [0.0] || return Float64[]
    else
      return [-q₀ / q₁]
    end
  end

  # Case where q(x) is indeed quadratic.
  rhs = sqrt(eps(Float64)) * q₁ * q₁
  if abs(q₀ * q₂) > rhs
    ρ = q₁ * q₁ - 4.0 * q₂ * q₀
    ρ < 0.0 && return Float64[]
    d = -0.5 * (q₁ + copysign(sqrt(ρ), q₁))
    roots = [d / q₂, q₀ / d]
  else
    # Ill-conditioned quadratic.
    roots = [-q₁ / q₂, 0.0]
  end

  # Perform a few Newton iterations to improve accuracy.
  for k = 1 : 2
    root = roots[k]
    for it = 1 : nitref
      q = (q₂ * root + q₁) * root + q₀
      dq = 2.0 * q₂ * root + q₁
      dq == 0.0 && continue
      root = root - q / dq
    end
    roots[k] = root
  end
  return roots
end


"""Given a trust-region radius `radius`, a vector `x` lying inside the
trust-region and a direction `d`, return `σ` > 0 such that

    ‖x + σ d‖ = radius

in the Euclidean norm. If known, ‖x‖² may be supplied in `xNorm2`.

If `flip` is set to `true`, `σ` > 0 is computed such that

    ‖x - σ d‖ = radius.
"""
function to_boundary(x :: Vector{Float64}, d :: Vector{Float64},
                     radius :: Float64; flip :: Bool=false, xNorm2 :: Float64=0.0)
  radius > 0 || error("radius must be positive")

  # σ is the positive root of the quadratic
  # ‖d‖² σ² + 2 xᵀd σ + (‖x‖² - radius²).
  xd = dot(x, d)
  flip && (xd = -xd)
  dNorm2 = dot(d, d)
  xNorm2 == 0.0 && (xNorm2 = dot(x, x))
  (xNorm2 <= radius * radius) || error(@sprintf("outside of the trust region: ‖x‖²=%7.1e, Δ²=%7.1e", xNorm2, radius * radius))
  roots = roots_quadratic(dNorm2, 2 * xd, xNorm2 - radius * radius)
  return maximum(roots)
end



"""Given a regularization parameter `regulα`, a vector `x` and a direction `d`, 
   return `σ` > 0 such that

    h(σ) ~ argmin_{α>0} h(α) = q(x + α d) +  ‖x + αd‖³/3regulα

in the Euclidean norm.
"""
function to_minimum(A :: LinearOperator, b :: Array{Float64,1},
                               x :: Vector{Float64}, 
                               p :: Vector{Float64}, Ap :: Vector{Float64}, pAp :: Float64,
                               α :: Float64,
                               regulα :: Float64)
    regulα > 0 || error("regulα must be positive")
    #println("regulα = ",regulα)

    Ax = copy(A * x)
    xAx = dot(x, Ax)
    xAp = dot(x, Ap)
    bx = dot(b, x)
    bp = dot(b, p)
    h = α -> 0.5*(xAx + 2.0*xAp*α + α^2*pAp) - bx - α*bp + norm(x+α*p)^3/(3*regulα)
    h0 = h(0)

    #dh = α -> xAp + α*pAp - bp + norm(x+α*p)*(x+α*p)'*p
    xx=dot(x,x)
    xp=dot(x,p)
    pp=dot(p,p)
    dh = α -> xAp + α*pAp - bp + sqrt(xx + 2*α*xp + α^2*pp)*(xp + α*pp)/regulα


#println("h0 = ",h(0),"  dh0 = ",dh(0))
    # search for a positive dh
    α = regulα
    while dh(α)<0
        α *= 5
    end
#println("h(α) = ",h(α),"  dh(α) = ",dh(α), "  α = ",α)

    a=0.0
    σ = α/2
assert((dh(a)<0) & (dh(α)>0))
    while (((α-a)/max(α,a))>1e-8) & (abs(dh(σ))>1e-8)
        if dh(σ)<0  a=σ
        else
            α = σ
        end
        σ = (a+α)/2
#println("  dh(σ) = ",dh(σ),"  dh(α) = ",dh(α),"  dh(a) = ",dh(a), "  σ = ",σ, "  α = ",α, "  a = ",a)
    end
#println("h(σ) = ",h(σ),"  dh(σ) = ",dh(σ), "  σ = ",σ)

    return σ
end
