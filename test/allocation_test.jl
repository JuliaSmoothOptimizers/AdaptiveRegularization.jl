using AdaptiveRegularization, LinearAlgebra, Test
using ADNLPModels, NLPModels, OptimizationProblems, Stopping

nlp = OptimizationProblems.ADNLPProblems.arglina()
n = nlp.meta.nvar
x = nlp.meta.x0
stp = NLPStopping(nlp)
H = hess(nlp, x)
g = grad(nlp, x)
ng = norm(g)
calls, max_calls = 0, 1000000

function alloc_preprocess(XData, H, g, ng, calls, max_calls, α)
  AdaptiveRegularization.preprocess!(XData, H, g, ng, calls, max_calls, α)
  return nothing
end

function alloc_solve_model(XData, H, g, ng, calls, max_calls, α)
  AdaptiveRegularization.solve_model!(XData, H, g, ng, calls, max_calls, α)
  return nothing
end

@testset "Allocation test in preprocess and solvemodel" begin
  for Data in (PDataKARC, PDataTRK, PDataST)
    XData = Data(S, T, n)
    @testset "Allocation test in preprocess with $(Data)" begin
      alloc_preprocess(XData, H, g, ng, calls, max_calls, α)
      al = @allocated alloc_preprocess(XData, H, g, ng, calls, max_calls, α)
      @test al == 0
    end

    @testset "Allocation test in with $(Data)" begin
      alloc_solve_model(XData, H, g, ng, calls, max_calls, α)
      al = @allocated alloc_solve_model(XData, H, g, ng, calls, max_calls, α)
      @test al == 0
    end
  end
end
