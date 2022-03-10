# using Compat
using LinearOperators


"Abstract  for statistics returned by a solver"
abstract type KrylovStats end

"Type for statistics returned by non-Lanczos solvers"
mutable struct SimpleStats <: KrylovStats
  solved :: Bool
  inconsistent :: Bool
  residuals :: Array{Float64,1}
  Aresiduals :: Array{Float64,1}
  status :: String
end

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
#
# First reformulate the TR case using the characteristic function to confirm it is
# equivalent to the usual implementation. A must be symmetric to define a quadratic objective
# q(x) = 0.5*x'*A*x - b'*x
#
#   JPD july 21 2016, Sherbrooke

include("krylov_aux.jl")

export cgARC

"""The conjugate gradient method to solve the symmetric linear system Ax=b.

The method does _not_ abort if A is not definite.
"""
function cgARC(A, b::Array{T, 1}; atol::Float64 = 1e-08, rtol::Float64 = 1e-06, itmax::Int = 0, regulα::Float64 = 1, verbose::Bool = false) where T <: Real
  n = size(b, 1);
  (size(A, 1) == n & size(A, 2) == n) || error("Inconsistent problem size");
  #isequal(triu(A)',tril(A)) || error("Must supply Hermitian matrix")
  regulα > 0.0  ||  error("Regularization must be strictly positive")

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
  verbose && @printf("%5d  %8.1e ", iter, rNorm)

  solved = rNorm <= ε;
  tired = iter >= itmax;
  on_boundary = false;
  status = "unknown";

  q = s ->  0.5 * dot(s, copy(A * s)) - dot(b, s)

  m = s ->  q(s) + norm(s)^3/(3*regulα)
  #hO = α -> m(x+α*p)


  while ! (solved || tired)
    Ap = copy(A * p);  # Bug in LinearOperators? A side effect spoils the computation without the copy.
    pAp = BLAS.dot(n, p, 1, Ap, 1);

    α = γ / pAp;

    # Compute step size to the min of the regularized model.
    σ = to_minimum(A, b, x, p, Ap, pAp, α, regulα)

    verbose && @printf("%8.1e  %7.1e  %7.1e\n", pAp, α, σ);

    # Move along p from x to the min if either
    # the next step leads farther than the min of the regularized model or
    # we have nonpositive curvature.
    if ((pAp <= 0.0) | (α > (σ/2)))
      α = σ
      on_boundary = true
      verbose &&println("at minimum after ",iter," CG iterations.")
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

  status = on_boundary ? "at min of the cubic regularization" : (tired ? "maximum number of iterations exceeded" : "solution good enough given atol and rtol")
  stats = SimpleStats(solved, false, rNorms, [], status);
  return (x, stats);
end
