module AdaptiveRegularization

# stdlib
using LinearAlgebra, SparseArrays
# JSO
using Krylov, LinearOperators, NLPModels, SparseMatricesCOO, SolverCore, SolverTools
# Stopping
using Stopping, StoppingInterface

# Selective includes.
include("./utils/hessian_rep.jl")
include("./utils/pdata_struct.jl")
include("./utils/increase_decrease.jl")

path = joinpath(dirname(@__FILE__), "SolveModel")
files = filter(x -> x[(end - 2):end] == ".jl", readdir(path))
for file in files
  include("SolveModel/" * file)
end

path = joinpath(dirname(@__FILE__), "PreProcess")
files = filter(x -> x[(end - 2):end] == ".jl", readdir(path))
for file in files
  include("PreProcess/" * file)
end

include("trarc_solver.jl") # structure for the solve!
include("main.jl") # main algo
include("solvers.jl")

export ALL_solvers, NLS_solvers

ALL_solvers = keys(solvers_const)
NLS_solvers = keys(solvers_nls_const)

export TRARC

"""
    TRARC(nlp; kwargs...)

Compute a local minimum of an unconstrained optimization problem using
trust-region (TR)/adaptive regularization with cubics (ARC) methods.

Some variants of TRARC are already implemented and listed in `AdaptiveRegularization.ALL_solvers`.

# Arguments
- `nlp::AbstractNLPModel`: the model solved, see `NLPModels.jl`.

The keyword arguments include
- `TR::ARTrustRegion`: structure with trust-region/ARC parameters, see [`SolverTools.jl`](https://github.com/JuliaSmoothOptimizers/SolverTools.jl). Default: `ARTrustRegion(T(10.0))`.
- `hess_type::Type{Hess}`: Structure used to handle the hessian. The possible values are: `HessDense`, `HessSparse`, `HessSparseCOO`, `HessOp`. Default: `HessOp`.
- `pdata_type::Type{ParamData}` Structure used for the preprocessing step. Default: `PDataKARC`.
- `robust::Bool`: `true` implements a robust evaluation of the model. Default: `true`.
- `verbose::Bool`: `true` prints iteration information. Default: `false`.
Additional `kwargs` are used for stopping criterion, see `Stopping.jl`.

# Output
The returned value is a `GenericExecutionStats`, see `SolverCore.jl`.

# Callback
The callback is called at each iteration.
The expected signature of the callback is `callback(nlp, solver, stats)`, and its output is ignored.
Changing any of the input arguments will affect the subsequent iterations.
In particular, setting `stats.status = :user` will stop the algorithm.
All relevant information should be available in `nlp` and `solver`.
Notably, you can access, and modify, the following:
- `solver.stp`: stopping object used for the algorithm;
- `solver.workspace`: additional allocations;
- `stats`: structure holding the output of the algorithm (`GenericExecutionStats`), which contains, among other things:
  - `stats.dual_feas`: norm of current gradient;
  - `stats.iter`: current iteration counter;
  - `stats.objective`: current objective function value;
  - `stats.status`: current status of the algorithm. Should be `:unknown` unless the algorithm has attained a stopping criterion. Changing this to anything will stop the algorithm, but you should use `:user` to properly indicate the intention.
  - `stats.elapsed_time`: elapsed time in seconds.


This implementation uses `Stopping.jl`. Therefore, it is also possible to used
        
    TRARC(stp; kwargs...)

which returns the `stp::NLPStopping` updated.

For advanced usage, first define a [`TRARCSolver`](@ref) to preallocate the
 memory used in the algorithm, and then call `solve!`:

    stats = solve!(solver, nlp)
    stats = solve!(solver, nlp, stats)

To choose a particular variant, the keyword arguments `hess_type` and `pdata_type` can be used as follows:
   `TRARCSolver(nlp; hess_type = ht, pdata_type = PDataKARC)`
the former specifying how to handle the Hessian information, i.e. matrix-free or sparse matrix, while
the latter specifies the type of adaptive approach, e.g. trust-region or ARC.

# References
This method unifies the implementation of trust-region and adaptive regularization with cubics as described in

    Dussault, J.-P. (2020).
    A unified efficient implementation of trust-region type algorithms for unconstrained optimization.
    INFOR: Information Systems and Operational Research, 58(2), 290-309.
    10.1080/03155986.2019.1624490

# Examples

```jldoctest
using AdaptiveRegularization, ADNLPModels
nlp = ADNLPModel(x -> 100 * (x[2] - x[1]^2)^2 + (x[1] - 1)^2, [-1.2; 1.0]);
stats = TRARC(nlp)

# output

"Execution stats: first-order stationary"
```

```jldoctest; output = false
using AdaptiveRegularization, ADNLPModels, SolverCore
nlp = ADNLPModel(x -> 100 * (x[2] - x[1]^2)^2 + (x[1] - 1)^2, [-1.2; 1.0]);
solver = TRARCSolver(nlp);
stats = solve!(solver, nlp)

# output

"Execution stats: first-order stationary"
```

```jldoctest; output = false
using AdaptiveRegularization, ADNLPModels, SolverCore
nlp = ADNLPModel(x -> 100 * (x[2] - x[1]^2)^2 + (x[1] - 1)^2, [-1.2; 1.0]);
solver = TRARCSolver(nlp);
stats = GenericExecutionStats(nlp)
stats = solve!(solver, nlp, stats)

# output

"Execution stats: first-order stationary"
```
"""
function TRARC end

function TRARC(nlp::AbstractNLPModel{T, S}; kwargs...) where {T, S}
  nlpstop = NLPStopping(nlp; optimality_check = (pb, state) -> norm(state.gx), kwargs...)
  stats = TRARC(nlpstop; kwargs...)
  return stats
end

function TRARC(nlp_stop::NLPStopping; kwargs...)
  nlp = nlp_stop.pb
  solver = TRARCSolver(nlp; kwargs...)
  stats = GenericExecutionStats(nlp)
  SolverCore.solve!(solver, nlp_stop, stats; kwargs...)
  return stats
end

for fun in union(keys(solvers_const), keys(solvers_nls_const))
  ht, pt, ka = merge(solvers_const, solvers_nls_const)[fun]
  @eval begin
    @doc (@doc TRARC) function $fun(nlpstop::NLPStopping; kwargs...)
      kw_list = Dict{Symbol, Any}()
      if $ka != ()
        for t in $ka
          push!(kw_list, t)
        end
      end
      merge!(kw_list, Dict(kwargs))
      TRARC(nlpstop; hess_type = $ht, pdata_type = $pt, kw_list...)
    end
  end
  @eval begin
    function $fun(nlp::AbstractNLPModel{T, S}; kwargs...) where {T, S}
      nlpstop = NLPStopping(nlp; optimality_check = (pb, state) -> norm(state.gx), kwargs...)
      stats = $fun(nlpstop; kwargs...)
      return stats
    end
  end
  @eval export $fun
end

end # module
