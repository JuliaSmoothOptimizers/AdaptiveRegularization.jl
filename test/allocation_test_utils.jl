using AdaptiveRegularization, LinearAlgebra, Test

T = Float64
S = Vector{T}
n = 1000

TR = TrustRegion(T(10))
α = T(100)

# Test increase / decrease

function alloc_decrease(XData, α, TR)
  AdaptiveRegularization.decrease(XData, α, TR)
  return nothing
end

function alloc_increase(XData, α, TR)
  AdaptiveRegularization.increase(XData, α, TR)
  return nothing
end

@testset "Allocation test in increase/decrease" begin
  for XData in (PDataKARC(S, T, n), PDataTRK(S, T, n), PDataST(S, T, n))
    @testset "Allocation test in decrease for $(typeof(XData))" begin
      alloc_decrease(XData, α, TR)
      al = @allocated alloc_decrease(XData, α, TR)
      @test al == 0
    end
    @testset "Allocation test in increase for $(typeof(XData))" begin
      alloc_increase(XData, α, TR)
      al = @allocated alloc_increase(XData, α, TR)
      @test al == 0
    end
  end
end

# Test hessian evaluation in-place

using NLPModelsTest
nlp = NLPModelsTest.BROWNDEN()
n = nlp.meta.nvar
x0 = nlp.meta.x0

function alloc_hessian(who, nlp, x0)
  AdaptiveRegularization.hessian!(who, nlp, x0)
  return nothing
end

@testset "Test in-place hessian allocations" begin
  @testset "$Workspace" for (Workspace, limit) in (
    (HessDense, 1952),
    (HessSparse, 944),
    (HessSparseCOO, 0),
    (HessOp, 960),
  )
    who = Workspace(nlp, n)
    alloc_hessian(who, nlp, x0)
    al = @allocated alloc_hessian(who, nlp, x0)
    @test al <= limit
    println("Allocations for $Workspace is $al.")
  end
end
