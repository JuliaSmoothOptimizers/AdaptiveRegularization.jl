using NLPModelsTest, Stopping, AdaptiveRegularization

nlp = BROWNDEN()
n, x0 = nlp.meta.nvar, nlp.meta.x0

TRnothing = TrustRegion(10.0, max_unsuccinarow = 2)

S, T = typeof(x0), eltype(x0)
PData = PDataKARC(S, T, n)

function alloc_AdaptiveRegularization(stp, PData, workspace, TRnothing)
  TRARC(stp, PData, workspace, TRnothing, solve_model = AdaptiveRegularization.solve_modelKARC)
  return nothing
end

@testset "Test TRARC allocations using $Hess" for (Hess, limit) in (
  (HessOp, 1904),
  (HessSparseCOO, 944),
)
  stp = NLPStopping(nlp)
  stp.meta.max_iter = 5
  workspace = AdaptiveRegularization.TRARCWorkspace(nlp, Hess)
  alloc_AdaptiveRegularization(stp, PData, workspace, TRnothing)
  al = @allocated alloc_AdaptiveRegularization(stp, PData, workspace, TRnothing)
  println("Allocations of TRARC using $Hess is $al .")
  @show al
end
