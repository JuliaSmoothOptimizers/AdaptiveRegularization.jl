# AdaptiveRegularization

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaSmoothOptimizers.github.io/AdaptiveRegularization.jl/stable)
[![In development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSmoothOptimizers.github.io/AdaptiveRegularization.jl/dev)
[![Build Status](https://github.com/JuliaSmoothOptimizers/AdaptiveRegularization.jl/workflows/Test/badge.svg)](https://github.com/JuliaSmoothOptimizers/AdaptiveRegularization.jl/actions)
[![Test workflow status](https://github.com/JuliaSmoothOptimizers/AdaptiveRegularization.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/JuliaSmoothOptimizers/AdaptiveRegularization.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Lint workflow Status](https://github.com/JuliaSmoothOptimizers/AdaptiveRegularization.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/JuliaSmoothOptimizers/AdaptiveRegularization.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/JuliaSmoothOptimizers/AdaptiveRegularization.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/JuliaSmoothOptimizers/AdaptiveRegularization.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaSmoothOptimizers/AdaptiveRegularization.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSmoothOptimizers/AdaptiveRegularization.jl)
[![DOI](https://zenodo.org/badge/DOI/FIXME)](https://doi.org/FIXME)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)
[![All Contributors](https://img.shields.io/github/all-contributors/JuliaSmoothOptimizers/AdaptiveRegularization.jl?labelColor=5e1ec7&color=c0ffee&style=flat-square)](#contributors)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)

AdaptiveRegularization is a solver for unconstrained nonlinear problems,

    min f(x)

It uses other [JuliaSmoothOptimizers](https://juliasmoothoptimizers.github.io/) packages for development.
In particular, [NLPModels.jl](https://github.com/JuliaSmoothOptimizers/NLPModels.jl) is used for defining the problem, and [SolverCore.jl](https://github.com/JuliaSmoothOptimizers/SolverCore.jl) for the output.

This package uses [`Stopping.jl`](https://github.com/SolverStoppingJulia/Stopping.jl) via `NLPStopping` to handle its workflow, you can also see [tutorials with `Stopping`](https://solverstoppingjulia.github.io/StoppingTutorials.jl) to learn more.

## Algorithm

The initial implementation of this package follows (Dussault, J.-P. 2020):

*Adaptive cubic regularization (ARC) and trust-region (TR) methods use modified linear systems to compute their steps. The modified systems consist in adding some multiple of the identity matrix (or a well-chosen positive definite matrix) to the Hessian to obtain a sufficiently positive definite linear system, the so called shifted system. This type of system was first proposed by Levenberg and Marquardt. Some trial and error is often involved to obtain a specified value for this shift parameter. We provide an efficient unified implementation to track the shift parameter; our implementation encompasses many ARC and TR variants.*

## References

> Dussault, J.-P. (2020).
> A unified efficient implementation of trust-region type algorithms for unconstrained optimization.
> INFOR: Information Systems and Operational Research, 58(2), 290-309.
> [10.1080/03155986.2019.1624490](https://doi.org/10.1080/03155986.2019.1624490)

> Dussault, J.-P., Migot, T. & Orban, D. (2023).
> Scalable adaptive cubic regularization methods.
> Mathematical Programming.
> [10.1007/s10107-023-02007-6](https://doi.org/10.1007/s10107-023-02007-6)

## How to Cite

If you use AdaptiveRegularization.jl in your work, please cite using the reference given in [CITATION.cff](https://github.com/JuliaSmoothOptimizers/AdaptiveRegularization.jl/blob/main/CITATION.cff).

## Contributing

If you want to make contributions of any kind, please first that a look into our [contributing guide directly on GitHub](docs/src/90-contributing.md) or the [contributing page on the website](https://JuliaSmoothOptimizers.github.io/AdaptiveRegularization.jl/dev/90-contributing/)

---

### Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->
