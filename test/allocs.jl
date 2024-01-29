"""
    @wrappedallocs(expr)
Given an expression, this macro wraps that expression inside a new function
which will evaluate that expression and measure the amount of memory allocated
by the expression. Wrapping the expression in a new function allows for more
accurate memory allocation detection when using global variables (e.g. when
at the REPL).
For example, `@wrappedallocs(x + y)` produces:
```julia
function g(x1, x2)
    @allocated x1 + x2
end
g(x, y)
```
You can use this macro in a unit test to verify that a function does not
allocate:
```
@test @wrappedallocs(x + y) == 0
```
"""
macro wrappedallocs(expr)
  argnames = [gensym() for a in expr.args]
  quote
    function g($(argnames...))
      @allocated $(Expr(expr.head, argnames...))
    end
    $(Expr(:call, :g, [esc(a) for a in expr.args]...))
  end
end

using AdaptiveRegularization
using NLPModelsTest, NLPModels, SolverCore

  @testset "Allocation tests" begin
    @testset "$symsolver with $ht" for symsolver in (:TRARCSolver,), ht in (HessOp, HessSparseCOO)
      for model in NLPModelsTest.nlp_problems
        nlp = eval(Meta.parse(model))()
        if unconstrained(nlp)
          solver = eval(symsolver)(nlp; hess_type = ht)
          stats = GenericExecutionStats(nlp)
          SolverCore.solve!(solver, nlp, stats)
          reset!(solver)
          reset!(nlp)
          al = @wrappedallocs SolverCore.solve!(solver, nlp, stats)
          @test al == 0
        end
      end
    end
  end

#=
using AdaptiveRegularization
using NLPModelsTest, NLPModels, SolverCore
problems = NLPModelsTest.nlp_problems
nlp = eval(Meta.parse(first(problems)))()
solver = TRARCSolver(nlp)
stats = GenericExecutionStats(nlp)

stp = solver.stp
PData = solver.meta
workspace = solver.workspace
TR = solver.TR
∇f = rand(nlp.meta.nvar)
norm_∇f = norm(∇f)
α = 1e-7
@code_warntype AdaptiveRegularization.get_hess(workspace.Hstruct)
@code_warntype AdaptiveRegularization.preprocess!(stp, PData, workspace, ∇f, norm_∇f, α)

using Profile, PProf
Profile.Allocs.clear()
NLPModels.reset!(nlp)
reset!(solver)
@time solve!(solver, nlp, stats)
NLPModels.reset!(nlp)
reset!(solver)
@time Profile.Allocs.@profile sample_rate=1 solve!(solver, nlp, stats)
PProf.Allocs.pprof(from_c = false)

=#