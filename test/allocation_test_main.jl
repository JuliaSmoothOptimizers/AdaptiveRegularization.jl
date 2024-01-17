using NLPModelsTest, Stopping, AdaptiveRegularization

nlp = BROWNDEN()
n, x0 = nlp.meta.nvar, nlp.meta.x0

TRnothing = TrustRegion(10.0, max_unsuccinarow = 2)

S, T = typeof(x0), eltype(x0)
PData = PDataKARC(S, T, n)

function alloc_AdaptiveRegularization(stp, solver, stats)
  solve!(solver, stp, stats, solve_model = AdaptiveRegularization.solve_modelKARC)
  return nothing
end

@testset "Test TRARC allocations using $Hess" for (Hess, limit) in
                                                  ((HessOp, 1904), (HessSparseCOO, 944))
  stp = NLPStopping(nlp)
  stp.meta.max_iter = 5
  solver = TRARCSolver(nlp, hess_type = Hess)
  stats = GenericExecutionStats(nlp)
  alloc_AdaptiveRegularization(stp, solver, stats)
  al = @allocated alloc_AdaptiveRegularization(stp, solver, stats)
  println("Allocations of TRARC using $Hess is $al .")
  @show al
end
