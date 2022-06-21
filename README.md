# ARCTR : A unified efficient implementation of trust-region type algorithms for unconstrained optimization

[![docs-stable][docs-stable-img]][docs-stable-url] [![docs-dev][docs-dev-img]][docs-dev-url] [![build-ci][build-ci-img]][build-ci-url] [![codecov][codecov-img]][codecov-url] [![release][release-img]][release-url]

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://vepiteski.github.io/ARCTR.jl/stable
[docs-dev-img]: https://img.shields.io/badge/docs-dev-purple.svg
[docs-dev-url]: https://vepiteski.github.io/ARCTR.jl/dev
[build-ci-img]: https://github.com/vepiteski/ARCTR.jl/workflows/CI/badge.svg?branch=main
[build-ci-url]: https://github.com/vepiteski/ARCTR.jl/actions
[codecov-img]: https://codecov.io/gh/vepiteski/ARCTR.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/vepiteski/ARCTR.jl
[release-img]: https://img.shields.io/github/v/release/vepiteski/ARCTR.jl.svg?style=flat-square
[release-url]: https://github.com/vepiteski/ARCTR.jl/releases

ARCTR is a solver for unconstrained nonlinear problems,

    min f(x)

It uses other [JuliaSmoothOptimizers](https://jso-docs.github.io) packages for development.
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

## How to Cite

If you use ARCTR.jl in your work, please cite using the format given in [CITATION.cff](https://github.com/vepiteski/ARCTR.jl/blob/main/CITATION.cff).  <!--https://citation-file-format.github.io/cff-initializer-javascript/#/ -->

## Installation

`pkg> add https://github.com/vepiteski/ARCTR.jl`

## Example

```julia
using ARCTR, ADNLPModels

# Rosenbrock
nlp = ADNLPModel(x -> 100 * (x[2] - x[1]^2)^2 + (x[1] - 1)^2, [-1.2; 1.0])
stats = ARCqKOp(nlp)
```

# Bug reports and discussions

If you think you found a bug, feel free to open an [issue](https://github.com/vepiteski/ARCTR.jl/issues).
Focused suggestions and requests can also be opened as issues. Before opening a pull request, start an issue or a discussion on the topic, please.

If you want to ask a question not suited for a bug report, feel free to start a discussion [here](https://github.com/JuliaSmoothOptimizers/Organization/discussions). This forum is for general discussion about this repository and the [JuliaSmoothOptimizers](https://github.com/JuliaSmoothOptimizers), so questions about any of our packages are welcome.
