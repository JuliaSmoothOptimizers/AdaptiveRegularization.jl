# Tutorial

We show here the basic features of the package.

```@example 1
using AdaptiveRegularization, ADNLPModels

# Rosenbrock
nlp = ADNLPModel(x -> 100 * (x[2] - x[1]^2)^2 + (x[1] - 1)^2, [-1.2; 1.0])
stats = ARCqKOp(nlp, verbose = true)
```

It is possible to access the number of evaluations of each function of the NLPModel API using the following:

```@example ex1
nobj = neval_obj(nlp) # return number of f call
ngra = neval_grad(nlp) # return number of gradient call
nhes = neval_hess(nlp) # return number of Hessian call
nhpr = neval_hprod(nlp) # return number of hessian-vector products

(nobj, ngra, nhes, nhpr)
```

These functions come from the NLPModel API defined in [NLPModels.jl](https://juliasmoothoptimizers.github.io/NLPModels.jl/dev/).
If you want to reset the internal counter, you just do `reset!(nlp)`.
