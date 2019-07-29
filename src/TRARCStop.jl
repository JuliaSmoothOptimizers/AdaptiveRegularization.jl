export TRARC2

function TRARC2(nlp 		:: AbstractNLPModel,
                nlp_stop 	:: NLPStopping;
                TR 			:: TrustRegion = TrustRegion(eltype(nlp.meta.x0)(10.0)),
			    c 			:: Combi = Combi{eltype(nlp.meta.x0)}(hessian_dense, PDataLDLt{eltype(nlp.meta.x0)}, solve_modelTRDiag, preprocessLDLt, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
				correction  :: Symbol = :None,
			    robust 		:: Bool = true,
                verbose 	:: Bool = false
                )

	T = eltype(nlp.meta.x0)
    # printstyled("On entre dans TRARC2 T = $T \n", color = :bold)

	nlp_at_x = nlp_stop.current_state
    hessian_rep, PData, solve_model, pre_process, decrease, params = extract(c)

    α = TR.α₀  # initial Trust Region size
	# @show typeof(α)
    xt, xtnext, d, Df = copy(nlp.meta.x0), copy(nlp.meta.x0), copy(nlp.meta.x0), 0.0
	# @show eltype(xt)
	xopt = xt
    λ = 1.0

    n = length(xt)
    ∇f = Array{T}(undef, n)
    ∇fnext = Array{T}(undef, n)

    ft = obj(nlp, xt)
	# @show typeof(ft)
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

	global xdemi = NaN * rand(n)

    while !OK
        PData = pre_process(nlp_at_x.Hx, ∇f, params, calls, nlp_stop.meta.max_eval)

		# printstyled("On a PData \n", color = :bold)

		# @show eltype(PData.L)
		# @show eltype(PData.D)

        if ~PData.OK
			@warn("Something wrong with PData")
			return nlp_at_x, nlp_stop.meta.optimal
		end

        success = false
		Ht = nothing
		# printstyled("On a succes = $success \n", color = :bold)
		# @show Ht

        while !success & !OK & (unsuccinarow < TR.max_unsuccinarow)
            try
                d, xdemi, λ = solve_model(nlp_stop, PData, α)
            catch
                println(" Problem in solve_model")
                return nlp_at_x, nlp_stop.meta.optimal
            end

            Δq = -(∇f + 0.5 * nlp_at_x.Hx * d)⋅d
			# @show eltype(Δq)

            if Δq < 0.0 println("*******   Ascent direction in SolveModel: Δq = $Δq")
                println("  g⋅d = $(∇f⋅d), 0.5 d'Hd = $(0.5*(nlp_at_x.Hx*d)⋅d)  α = $α  λ = $λ")
				#@bp
                #try println(" cond(H) = $(cond(full(H)))") catch println("sparse hessian, no cond") end
                #D, Q = eig(full(H))
                #lm = findmin(D); lM = findmax(D)
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
			# printstyled("On a slope = $slope \n", color = :bold)
            # xtnext = xt + d
			if !(true in isnan.(xdemi))
				xtnext = xdemi + d
			else
				xtnext = xt + d
			end

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
			# printstyled("On a r = $r \n", color = :bold)
			# printstyled("On a good_grad = $good_grad \n", color = :bold)
			# printstyled("On a ∇fnext = $∇fnext \n", color = :bold)


            if r < TR.acceptance_threshold
                verbose && display_failure(iter, ftnext, λ, α)
	        	unsucc = unsucc + 1
	        	unsuccinarow = unsuccinarow + 1
	        	α = decrease(PData, α, TR)
                # fbidon = obj(nlp, xt)
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
				verysucc += 1
                Ht = hessian_rep(nlp, xt)
		        if r > TR.increase_threshold
		            α = increase(PData, α, TR)
		            verbose && display_v_success(iter, ft, norm_∇f, λ, α)
		        else
	                if r < TR.reduce_threshold
	                    α = decrease(PData, α, TR)
	                end
		            verbose && display_success(iter, ft, norm_∇f, λ, α)
		            succ += 1
		        end
	    	end
        end # while !succes

		OK = update_and_stop!(nlp_stop, x = xt, fx = ft, gx = ∇f, Hx = Ht)
		nlp_stop.meta.nb_of_stop = iter
		# @show eltype(nlp_at_x.x)
		calls = [nlp.counters.neval_obj,  nlp.counters.neval_grad, nlp.counters.neval_hess, nlp.counters.neval_hprod]
    end # while !OK

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
