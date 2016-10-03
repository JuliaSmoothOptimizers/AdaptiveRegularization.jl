# A standard implementation of the Conjugate Gradient method.
# The only non-standard point about it is that it does not check
# that the operator is definite.
# It is possible to check that the system is inconsistent by
# monitoring ‖p‖, which would cost an extra norm computation per
# iteration.
#
# Dominique Orban, <dominique.orban@gerad.ca>
# Salt Lake City, UT, March 2015.
#
# Try to adapt a stopping criterion based on regulα for the ARCq.
# First reformulate the TR case using the characteristic function to confirm it is 
# equivalent to the usual implementation. A must be symmetric to define a quadratic objective
# q(x) = 0.5*x'*A*x - b'*x
#
#   JPD july 21 2016, Sherbrooke

include("krylov_aux.jl")

export cgTR

# Methods for various argument types.
#include("cg_methods.jl")

"""The conjugate gradient method to solve the symmetric linear system Ax=b.

The method does _not_ abort if A is not definite.
"""
function cgTR{T <: Real}(A :: LinearOperator, b :: Array{T,1};
                       atol :: Float64=1.0e-8, rtol :: Float64=1.0e-6, itmax :: Int=0,
                       radius :: Float64=0.0, verbose :: Bool=false)

  n = size(b, 1);
  (size(A, 1) == n & size(A, 2) == n) || error("Inconsistent problem size");
  verbose && @printf("CG: system of %d equations in %d variables\n", n, n);

  # Initial state.
  x = zeros(n)
  x̂ = copy(x)

  γ = dot(b, b);
  γ == 0 && return x;
  r = copy(b);
  p = copy(r);

  σ = 0.0

  iter = 0;
  itmax == 0 && (itmax = 2 * n);

  rNorm = sqrt(γ);
  rNorms = [rNorm;];
  ε = atol + rtol * rNorm;
  verbose && @printf("%5d  %8.1e ", iter, rNorm);

  solved = rNorm <= ε;
  tired = iter >= itmax;
  on_boundary = false;
  status = "unknown";

  q = s ->  0.5 * dot(s, copy(A * s)) - dot(b, s) 
  m = s -> norm(s) > radius ? Inf : q(s)

  while ! (solved || tired)
    Ap = copy(A * p);  # Bug in LinearOperators? A side effect spoils the computation without the copy.
    pAp = BLAS.dot(n, p, 1, Ap, 1);

    α = γ / pAp;

#    # Compute step size to boundary if applicable.
#    σo = radius > 0.0 ? to_boundary(x, p, radius) : α
#
#    verbose && @printf("%8.1e  %7.1e  %7.1e\n", pAp, α, σo);
#
#    # Move along p from x to the boundary if either
#    # the next step leads outside the trust region or
#    # we have nonpositive curvature.
#    if (radius > 0.0) & ((pAp <= 0.0) | (α > σo))
#      αo = σo
#      on_boundary = true
#    end
#    verbose && @printf("                      %8.1e  %7.1e  %7.1e\n", pAp, αo, σo);
          
    if pAp <= 0.0  α=Inf;  end # negative curvature

    x̂ = copy(x)
    BLAS.axpy!(n,  α,  p, 1, x̂, 1);  # Faster than x̂ = x̂ + α * p;
    if radius > 0.0   # check now for the boundary via the model function m
       if m(x̂) > 0.9*m(x)
#           # minimize m(x+σ*p) ; we know the optimum lies on the boundary
#           #σ = to_boundary(x, p, radius)
           α = to_boundary(x, p, radius)
           on_boundary = true
       end
    end
    verbose && @printf("    %8.1e  %7.1e  %7.1e\n", pAp, α, σ);
    BLAS.axpy!(n,  α,  p, 1, x, 1);  # Faster than x = x + σ * p;
    BLAS.axpy!(n, -α, Ap, 1, r, 1);  # Faster than r = r - α * Ap;
    γ_next = BLAS.dot(n, r, 1, r, 1);
    rNorm = sqrt(γ);
    push!(rNorms, rNorm);

    solved = (rNorm <= ε) | on_boundary;
    tired = iter >= itmax;

    if !solved
      β = γ_next / γ;
      γ = γ_next;
      BLAS.scal!(n, β, p, 1)
      BLAS.axpy!(n, 1.0, r, 1, p, 1);  # Faster than p = r + β * p;
    end
    iter = iter + 1;
    verbose && @printf("%5d  %8.1e ", iter, rNorm);
  end
  verbose && @printf("\n");

  status = on_boundary ? "on trust-region boundary" : (tired ? "maximum number of iterations exceeded" : "solution good enough given atol and rtol")
  stats = SimpleStats(solved, false, rNorms, [], status);
  return (x, stats);
end
