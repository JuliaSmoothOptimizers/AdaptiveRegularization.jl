# Scalable Adaptive Cubic Regularization Methods

This folder contains a set of scripts to reproduce the results in

> Dussault, J.-P., Migot, T. & Orban, D. (2023).
> Scalable adaptive cubic regularization methods.
> Mathematical Programming.
> [10.1007/s10107-023-02007-6](https://doi.org/10.1007/s10107-023-02007-6)

There are three scripts:

- `make_problem_cutest_list.jl`: generate a data file with the list of problems solved;
- `benchmark_arctr.jl`: run the benchmark comparing ST_TR and ARCqK, and store the result in a JLD2 file;
- `figure_solvers.jl`: generates a set of plots.

The result file `2022-05-16_ST_TROp_ARCqKOpShift05_cutest_277_1000000.jld2` can be accessed without reproducing the benchmark.
