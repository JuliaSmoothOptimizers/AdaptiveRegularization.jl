```@meta
CurrentModule = AdaptiveRegularization
```

# AdaptiveRegularization

Documentation for [AdaptiveRegularization](https://github.com/JuliaSmoothOptimizers/AdaptiveRegularization.jl).

AdaptiveRegularization is a solver for unconstrained nonlinear problems,

    min f(x)

It uses other [JuliaSmoothOptimizers](https://juliasmoothoptimizers.github.io/) packages for development.
In particular, [NLPModels.jl](https://github.com/JuliaSmoothOptimizers/NLPModels.jl) is used for defining the problem, and [SolverCore.jl](https://github.com/JuliaSmoothOptimizers/SolverCore.jl) for the output.

This package uses [`Stopping.jl`](https://github.com/SolverStoppingJulia/Stopping.jl) via `NLPStopping` to handle its workflow, you can also see [tutorials with `Stopping`](https://solverstoppingjulia.github.io/StoppingTutorials.jl) to learn more.

## Algorithm

The initial implementation of this package follows (Dussault, J.-P. 2020):

*Adaptive cubic regularization (ARC) and trust-region (TR) methods use modified linear systems to compute their steps. The modified systems consist in adding some multiple of the identity matrix (or a well-chosen positive definite matrix) to the Hessian to obtain a sufficiently positive definite linear system, the so called shifted system. This type of system was first proposed by Levenberg and Marquardt. Some trial and error is often involved to obtain a specified value for this shift parameter. We provide an efficient unified implementation to track the shift parameter; our implementation encompasses many ARC and TR variants.

## Contributors

```@raw html
<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
```
