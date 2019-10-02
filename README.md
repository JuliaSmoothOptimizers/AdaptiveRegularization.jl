# ARCTR
[![Build Status](https://travis-ci.org/Goysa2/ARCTR.jl.svg?branch=master)](https://travis-ci.org/Goysa2/ARCTR.jl)

[![Coverage Status](https://coveralls.io/repos/github/Goysa2/ARCTR.jl/badge.svg?branch=master)](https://coveralls.io/github/Goysa2/ARCTR.jl?branch=master)

Several ARC and TR optimization solvers.


## Purpose
This package implement several Trust-Region and ARC methods to solve the unconstrained problem
<p align="center">
<b><i> min f(x) </b></i>
</p>

where <b>f</b> is a twice continuously differentiabe function.

### Trust-Region subproblem
Multiple ways are provided to solve the Trust-Region subproblem, which is to find a direction d(λ):
<p align="center">
<b><i> (∇²f(x)+λI)d(λ) = -∇f(x) </b></i>
</p>

such that ||d(λ)|| ⩽ Δ, with Δ being the size of the trust region.



## Installing
The optimal use of this package is through the <b>State</b> and <b>Stopping</b> packages.
```
] add https://github.com/Goysa2/State.jl
] add https://github.com/Goysa2/Stopping.jl
```

Although it is possible to use ARCTR in a self contained manner, support and update for future version of Julia will be garanteed only for usage with State and Stopping.

```
] add https://github.com/Goysa2/ARCTR.jl
```

## Usage example
Let's solve a famous problem: minimize the Rosenbrock function.
```
function rosenbrock(x)
	n = 2; m = 2;
	f = []
	push!(f, 10 * (x[2]-x[1]^2))
	push!(f, (x[1]-1))
	return sum(f[i]^2 for i=1:m)
end
```

We have to use <b>NLPModels.jl</b> to put our problems in a structure our algorithms can understand. Our starting point will be <b>[-1.2, 1.0]</b>. We also have to create our State and Stopping objects for this problem.

```
nlp = ADNLPModel(rosenbrock, [-1.2, 1.0]);
nlpstop = NLPStopping(nlp, Stopping.unconstrained, NLPAtX([-1.2, 1.0]));
```

We offer multiple ways to now solve the problem:
```
final_state, optimal = TRLDLt(nlp, nlpstop, verbose = true)
final_state, optimal = TRLDLt_abs(nlp, nlpstop, verbose = true)
final_state, optimal = TRSpectral(nlp, nlpstop, verbose = true)
final_state, optimal = TRSpectral_abs(nlp, nlpstop, verbose = true)
final_state, optimal = ARCSpectral(nlp, nlpstop, verbose = true)
final_state, optimal = ARCLDLt(nlp, nlpstop, verbose = true)
final_state, optimal = ARCqKOp(nlp, nlpstop, verbose = true)
```

The `final_state` provide information at the last iteration and `optimal` is a boolean value saying if the problem has reached an optimal solution or not.

### High-order correction
A limited selection of our methods offer the option to add an high-order correction if we reduce to Newton method's. The high-order correction we offer are Shamanskii, Chebyshev, Halley and SuperHalley. More documentation will be provided when those methods are more developed. 


## Long term goals
	* Introduce dynamic precision in algorithms
	* Introduce high order correction in an efficient manner.
	* Provide more documentation and examples
	* Making the package more self contained.
