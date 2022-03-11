# test all solvers with the well known Woods test function
nlp = woods(n = 10)

solver = ALL_solvers[21]
(x, f, gNorm, iter, optimal, tired, status) = solver(nlp, verbose = true)

reset!(nlp);
