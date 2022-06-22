using ARCTR, LinearAlgebra, Test

T = Float64
S = Vector{T}
n = 1000

TR = TrustRegion(T(10))
α = T(100)

for XData in (
    PDataKARC(S, T, n),
    PDataTRK(S, T, n),
    PDataST(S, T, n),
    # PDataSpectral(S, T, n),
)
    @testset "Allocation test in ARCTR.decrease for $(typeof(XData))" begin
        alloc_decrease() = @allocated ARCTR.decrease(XData, α, TR)
        alloc_decrease()
        @test alloc_decrease() <= 16
    end
    @testset "Allocation test in ARCTR.decrease for $(typeof(XData))" begin
        alloc_increase() = @allocated ARCTR.increase(XData, α, TR)
        alloc_increase()
        @test (@allocated alloc_increase()) <= 16
    end
end

using ADNLPModels, NLPModels, OptimizationProblems, Stopping

nlp = OptimizationProblems.ADNLPProblems.arglina()
n = nlp.meta.nvar
x = nlp.meta.x0
stp = NLPStopping(nlp)
H = hess(nlp, x)
g = grad(nlp, x)
ng = norm(g)
calls, max_calls = 0, 1000000

for (Data, solve, limit_solve, limit_preprocess) in (
    #:solve_diag,
    #:solve_diagTR,
    # (PDataSpectral(S, T, n), :solve_modelARCDiag),
    # (PDataSpectral(S, T, n), :solve_modelARCDiagAbs),
    (PDataKARC, :solve_modelKARC, 96, 4488),
    (PDataTRK, :solve_modelTRK, 96, 2400),
    (PDataST, :solve_modelST_TR, 544, 0),
    #(PDataSpectral, :solve_modelTRDiag, 5408, 279328),
    #(PDataSpectral, :solve_modelTRDiagAbs, 5408, 279328),
)
    @testset "Allocation test in preprocess with $(Data)" begin
        XData = Data(S, T, n)
        alloc_preprocess(XData, H, g, ng, calls, max_calls, α) =
            @allocated ARCTR.preprocess(XData, H, g, ng, calls, max_calls, α)
        alloc_preprocess(XData, H, g, ng, calls, max_calls, α)
        @test alloc_preprocess(XData, H, g, ng, calls, max_calls, α) <= limit_preprocess
        @show alloc_preprocess(XData, H, g, ng, calls, max_calls, α)
    end

    @testset "Allocation test in $solve with $(Data)" begin
        XData = Data(S, T, n)
        alloc_solve_model(H, g, ng, stp, XData, α) =
            @allocated ARCTR.eval(solve)(H, g, ng, stp, XData, α)
        alloc_solve_model(H, g, ng, stp, XData, α)
        @test alloc_solve_model(H, g, ng, stp, XData, α) <= limit_solve
        @show alloc_solve_model(H, g, ng, stp, XData, α)
    end
end
