# test all solvers with the well known Woods test function
nlp = woods(n = 10)

solver = ALL_solvers[25]
(x, f, gNorm, iter, optimal, tired, status) = solver(nlp, verbose = false);
@printf(
    "%-15s  %8d  %9.2e  %7.1e  %5d  %5d  %6d  %s\n",
    nlp.meta.name,
    nlp.meta.nvar,
    f,
    gNorm,
    nlp.counters.neval_obj,
    nlp.counters.neval_grad,
    nlp.counters.neval_hprod,
    status
)    #stats = run_solver(solver, model, verbose=false)
#@test (all([stats...] .>= 0))
reset!(nlp)
solver = ALL_solvers[26]
(x, f, gNorm, iter, optimal, tired, status) = solver(nlp, verbose = false);
@printf(
    "%-15s  %8d  %9.2e  %7.1e  %5d  %5d  %6d  %s\n",
    nlp.meta.name,
    nlp.meta.nvar,
    f,
    gNorm,
    nlp.counters.neval_obj,
    nlp.counters.neval_grad,
    nlp.counters.neval_hprod,
    status
)    #stats = run_solver(solver, model, verbose=false)

#reset!(nlp);
