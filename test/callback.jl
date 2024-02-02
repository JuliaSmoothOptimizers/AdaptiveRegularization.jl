@testset "Test callback" begin
  f(x) = (x[1] - 1)^2 + 4 * (x[2] - x[1]^2)^2
  nlp = ADNLPModel(f, [-1.2; 1.0])
  function cb(nlp, solver, stats)
    if stats.iter == 8
      stats.status = :user
    end
  end
  stats = TRARC(nlp, callback = cb)
  @test stats.iter == 8
end

@testset "Test callback for NLS" begin
  F(x) = [x[1] - 1; 2 * (x[2] - x[1]^2)]
  nls = ADNLSModel(F, [-1.2; 1.0], 2)
  function cb(nlp, solver, stats)
    if stats.iter == 8
      stats.status = :user
    end
  end

  stats = TRARC(nls, callback = cb)
  @test stats.iter == 8
end
