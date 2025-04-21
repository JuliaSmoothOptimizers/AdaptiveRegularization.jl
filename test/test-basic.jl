@testset "Testing NLP solvers" begin
  @testset "$name" for name in ALL_solvers
    solver = eval(name)
    unconstrained_nlp(solver)
    multiprecision_nlp(solver, :unc, precisions = (Float32, Float64))
  end
end

@testset "Testing NLS solvers" begin
  @testset "$name" for name in union(ALL_solvers, NLS_solvers)
    solver = eval(name)
    unconstrained_nls(solver)
    multiprecision_nls(solver, :unc, precisions = (Float32, Float64))
  end
end

@testset "Invidivual tests on all solvers" begin
  global nbsolver = 0
  verbose = false # turn true for debug
  @testset "Testing $solver individually" for solver in ALL_solvers
    global nbsolver += 1
    nlp = extrosnb(n = 2)
    nlpstop = NLPStopping(nlp)
    verbose && println(nbsolver, "  ", solver)
    eval(solver)(nlpstop, verbose = verbose)
    final_nlp_at_x, optimal = nlpstop.current_state, nlpstop.meta.optimal
    @test optimal
    NLPModels.reset!(nlp)
    stats = eval(solver)(nlp, verbose = verbose)
    @test stats.status == :first_order
    NLPModels.reset!(nlp)
  end
end
