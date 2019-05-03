export TRARC2_Shamanskii

function TRARC2_Shamanskii(nlp 		:: AbstractNLPModel,
                		   nlp_stop :: NLPStopping;
                		   TR 		:: TrustRegion = TrustRegion(10.0),
			    		   c 		:: Combi =  Combi(hessian_dense, PDataLDLt, solve_modelTRDiag, preprocessLDLt, decreaseFact, Tparam()),
			    		   robust 	:: Bool = true,
                		   verbose 	:: Bool = true
                		   )

	nlp_at_x = nlp_stop.current_state

    hessian_rep, PData, solve_model, pre_process, decrease, params = extract(c)


    α = TR.α₀  # initial Trust Region size
    xt, xtnext, d, Df = copy(nlp.meta.x0), copy(nlp.meta.x0), copy(nlp.meta.x0), 0.0
    xopt = xt
    λ = 1.0

    n = length(xt)
    ∇f = Array{Float64}(undef, n)
    ∇fnext = Array{Float64}(undef, n)

    ft = obj(nlp, xt)

	OK = update_and_start!(nlp_stop, x = xt, fx = ft, gx = ∇f, g0 = ∇f)


    fopt = ft
    grad!(nlp, xt, ∇f)
	OK = update_and_start!(nlp_stop, x = xt, fx = ft, gx = ∇f, g0 = ∇f)

	norm_∇f = norm(nlp_at_x.gx)
    norm_∇f0 = norm_∇f
    ∇fopt = ∇f
    norm_∇fopt = norm_∇f
	!OK && update!(nlp_at_x, Hx = hessian_rep(nlp, xt))


    ftnext = ft
    iter = 0

    verbose && display_header_iterations()
    verbose && display_success(iter, ftnext, norm_∇f0, 0.0, α)

    succ, unsucc, verysucc, unsuccinarow = 0, 0, 0, 0

    max_unsuccinarow = TR.max_unsuccinarow

    calls = [0, 0, 0, 0]

    while !OK # ~finish

        PData = pre_process(nlp_at_x.Hx, ∇f, params, calls, nlp_stop.meta.max_eval)

        if ~PData.OK
			@warn("Something wrong with PData")
			return nlp_at_x, nlp_stop.meta.optimal
		end

        success = false

        while ~success & !OK & (unsuccinarow < TR.max_unsuccinarow)
            try
                d, λ = solve_model(PData, α)
            catch
                println(" Problem in solve_model")

                return nlp_at_x, nlp_stop.meta.optimal
            end

			printstyled("avant la modification \n", color = :cyan)
			@show d
			xtemp = xt + d
			gtemp = grad(nlp, xtemp)
			d = nlp_at_x.Hx \ gtemp
			printstyled("après la modification \n", color = :cyan)
			@show d

            Δq = -(∇f + 0.5 * nlp_at_x.Hx * d)⋅d
			# printstyled("le Δq original est $(-(∇f + 0.5 * nlp_at_x.Hx * d)⋅d) \n", color = :cyan)
			# Δq  = nlp_at_x.gx' * d # just a test

			# printstyled("autre test de vérification pour direction de monté \n", color = :cyan)
			# @show nlp_at_x.gx' * d

            if Δq < 0.0 println("*******   Ascent direction in SolveModel: Δq = $Δq")
                println("  g⋅d = $(∇f⋅d), 0.5 d'Hd = $(0.5*(nlp_at_x.Hx*d)⋅d)  α = $α  λ = $λ")
				#@bp
                #try println(" cond(H) = $(cond(full(H)))") catch println("sparse hessian, no cond") end
                #D,Q=eig(full(H))
                #lm=findmin(D);lM=findmax(D)
                #println("λ min (H) = $lm  λ max (H) = $lM")
                #try printtln(" cond(L) = $(cond(PData.L))") catch println("Spectral") end
                #try println(" ||H - reconstructH|| = $(norm(H-reconstructH(PData)))")
                #catch
                #    println(" reconstruct H no available ")
                #end
                println("  calls  $calls   iters $verysucc; $succ; $unsucc")

				return nlp_at_x, nlp_stop.meta.optimal
            end
            slope = ∇f ⋅ d
            xtnext = xt + d

            iter += 1

            try
                ftnext = obj(nlp, xtnext)
            catch
                ftnext = Inf;
            end
            if isnan(ftnext)
                ftnext = Inf
            end

            Δf = ft - ftnext

            r, good_grad, ∇fnext = compute_r(nlp, ft, Δf, Δq, slope, d, xtnext, ∇fnext, robust)

            if r < TR.acceptance_threshold
                verbose && display_failure(iter, ftnext, λ, α)
	        	unsucc = unsucc + 1
	        	unsuccinarow = unsuccinarow + 1
	        	α = decrease(PData, α, TR)
                fbidon = obj(nlp, xt)
	    	else
	        	success = true

	        	unsuccinarow = 0
	        	xt = copy(xtnext)
	        	ft = ftnext

	        	if good_grad
					∇f = copy(∇fnext)
            	else
					grad!(nlp, xt, ∇f);
                end

                norm_∇f = stop_norm(∇f)
                H = hessian_rep(nlp, xt)
		        if r > TR.increase_threshold
		            α = increase(PData, α, TR)
		            verbose && display_v_success(iter, ft, norm_∇f, λ, α)
	                    verysucc += 1
		        else
	                    if r < TR.reduce_threshold
	                        α = decrease(PData, α, TR)
	                    end
		            verbose && display_success(iter, ft, norm_∇f, λ, α)
		            succ += 1
		        end
	    	end
        end

		OK = update_and_stop!(nlp_stop, x = xt, fx = ft, gx = ∇f, Hx = hessian_rep(nlp, xt))
        calls = [nlp.counters.neval_obj,  nlp.counters.neval_grad, nlp.counters.neval_hess, nlp.counters.neval_hprod]
    end

    xopt = xt
    fopt = ft
    ∇fopt = ∇f
    norm_∇fopt = norm_∇f

    niter = verysucc + succ + unsucc;
    verbose && @printf(" %10i unsuccessful iterations \n", unsucc)
    verbose && @printf(" %10i successful iterations \n", succ)
    verbose && @printf(" %10i very successful iterations \n", verysucc)
    verbose && @printf(" %10i total successful iterations \n", verysucc + succ)
    verbose && @printf(" %10i total iterations \n", niter)

    calls = [nlp.counters.neval_obj, nlp.counters.neval_grad, nlp.counters.neval_hess, nlp.counters.neval_hprod]

    total_calls = nlp.counters.neval_obj + nlp.counters.neval_grad + n * nlp.counters.neval_hess + nlp.counters.neval_hprod
     verbose && @printf(" %10i total f %10i g %10i H calls  %10i Hv calls  \n", nlp.counters.neval_obj, nlp.counters.neval_grad, nlp.counters.neval_hess, nlp.counters.neval_hprod)


    if OK && verbose  @printf("\n     Final cost(x) = %f  ||g|| = %9.2e  nb of f/g/H evals = %i  converged\n", fopt, norm_∇fopt, total_calls)
    elseif verbose
        @printf("\n     Final cost(x) = %f  ||g|| = %9.2e  nb of f/g/H evals  = %i  NOT converged\n\n", fopt, norm_∇fopt, total_calls)
    end;

    return nlp_at_x, nlp_stop.meta.optimal

end
