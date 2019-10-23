export TRARC

function TRARC(nlp 		:: AbstractNLPModel,
                nlp_stop 	:: NLPStopping;
                TR 			:: TrustRegion = TrustRegion(eltype(nlp.meta.x0)(10.0)),
			    c 			:: Combi = Combi{eltype(nlp.meta.x0)}(hessian_dense, PDataLDLt{eltype(nlp.meta.x0)}, solve_modelTRDiag, preprocessLDLt, decreaseFact, Tparam{eltype(nlp.meta.x0)}()),
			    robust 		:: Bool = true,
                verbose 	:: Bool = false
                )

	T = eltype(nlp_stop.current_state.x)
	# @show nlp.meta.x0
	dTR = nothing
	dN = nothing

	nlp_at_x = nlp_stop.current_state
    hessian_rep, PData, solve_model, pre_process, decrease, params = extract(c)
	# println("on a extract")
	# @show hessian_rep
	# @show PData
	# @show solve_model
	# @show pre_process
	# @show decrease
	# @show params

    α = TR.α₀  # initial Trust Region size
    xt, xtnext, d, Df = copy(nlp_at_x.x), copy(nlp_at_x.x), copy(nlp_at_x.x), 0.0
	xopt = xt
    λ = 1.0

    n = length(xt)
    ∇f = Array{T}(undef, n)
    ∇fnext = Array{T}(undef, n)
	# ∇fnext

    ft = obj(nlp, xt)
    fopt = ft
    grad!(nlp, xt, ∇f)
	OK = update_and_start!(nlp_stop, x = xt, fx = ft, gx = ∇f, g0 = copy(∇f))

	norm_∇f = norm(nlp_at_x.gx)
    norm_∇f0 = norm_∇f
    ∇fopt = ∇f
    norm_∇fopt = norm_∇f
	!OK && State.update!(nlp_at_x, Hx = hessian_rep(nlp, xt))


	# global cgtol = 1.0

    ftnext = ft
    iter = 0

    # verbose && display_header_iterations()
    # verbose && display_success(iter, ftnext, norm_∇f0, 0.0, α)
	@info log_header([:iter, :f, :nrm_g, :λ,  :status, :α, :nrm_dtr, :nrm_dHO, :f_p_dTR, :f_p_dHO, :dir, :ΔqN], [Int64, T, T, T, String, T, T, T, T, T, String, T])
	@info log_row(Any[iter, nlp_at_x.fx, norm_∇f, 0.0, "First iteration", α])

	# @info log_header([:iter, :f, :nrm_g, :λ,  :status, :α, :nrm_dtr], [Int64, T, T, T, String, T, T])
	# @info log_row(Any[iter, nlp_at_x.fx, norm_∇f, 0.0, "First iteration", α])


    succ, unsucc, verysucc, unsuccinarow = 0, 0, 0, 0

    max_unsuccinarow = TR.max_unsuccinarow

    calls = [0, 0, 0, 0]

	# global xdemi = NaN * rand(n)
	global dTR = nothing
	global dHO = nothing
	global dir = nothing


    while !OK
        PData = pre_process(nlp_at_x.Hx, ∇f, params, calls, nlp_stop.meta.max_eval)

        if ~PData.OK
			@warn("Something wrong with PData")
			return nlp_at_x, nlp_stop.meta.optimal
		end

        success = false
		Ht = nothing

        while !success & !OK & (unsuccinarow < TR.max_unsuccinarow)
			# printstyled("yoy! \n", color = :green)
			try
				# @show solve_model
            	dTR, dHO, λ = solve_model(nlp_stop, PData, α)
            catch
            	println(" Problem in solve_model")
            	return nlp_at_x, nlp_stop.meta.optimal
            end
			# @show eltype(dTR)
			if dTR == dHO
				dir = "dTR == dHO"
			else
				dir = "dTR != dHO"
			end

			d = copy(dTR)
			# @show d[1]
			# @show d[2]

            Δq = -(∇f + 0.5 * nlp_at_x.Hx * d)⋅d

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
			xt_original = copy(xt)
			xtnext = xt .+ d
			# @show xtnext[1]
			# @show xtnext[2]

            iter += 1

            try
                ftnext = obj(nlp, xtnext)
            catch
                ftnext = Inf;
            end

			# printstyled("on a ftnext ↑ \n", color = :yellow)

            if isnan(ftnext)
                ftnext = Inf
            end

            Δf = ft - ftnext

            r, good_grad, ∇fnext = compute_r(nlp, ft, Δf, Δq, slope, d, xtnext, ∇fnext, robust)

            if r < TR.acceptance_threshold
				# printstyled("on est unsucessful \n", color = :yellow)
				@info log_row(Any[iter, nlp_at_x.fx, norm_∇f, λ, "U", α, norm(d), norm(dHO), obj(nlp, xt_original .+ d), obj(nlp, xt_original .+ dHO), dir, Δq])
				# @info log_row(Any[iter, nlp_at_x.fx, norm_∇f, λ, "Unsuccessful", α, norm(d)])
                # verbose && display_failure(iter, ftnext, λ, α)
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
					# printstyled("on est very sucessful \n", color = :yellow)
		            α = increase(PData, α, TR)
					# @show iter
					# @show ft
					# @show norm_∇f
					# @show λ
					# @show α
					@info log_row(Any[iter, ft, norm_∇f, λ, "V", α, norm(d), norm(dHO), obj(nlp, xt_original + dTR), obj(nlp, xt_original + dHO), dir, Δq])
					# @info log_row(Any[iter, ft, norm_∇f, λ, "Very Successful", α, norm(d)])
		            # verbose && display_v_success(iter, ft, norm_∇f, λ, α)
		        else
	                if r < TR.reduce_threshold
	                    α = decrease(PData, α, TR)
	                end
					# printstyled("on est sucessful \n", color = :yellow)
					@info log_row(Any[iter, ft, norm_∇f, λ, "S", α, norm(d), norm(dHO),  obj(nlp, xt_original + dTR), obj(nlp, xt_original + dHO), dir, Δq])
					# @info log_row(Any[iter, ft, norm_∇f, λ, "Successful", α, norm(d)])
		            # verbose && display_success(iter, ft, norm_∇f, λ, α)
		            succ += 1
		        end
	    	end
        end # while !succes

		# @info log_row([iter, nlp_at_x.fx, norm_∇f, λ, α])

		OK = update_and_stop!(nlp_stop, x = xt, fx = ft, gx = ∇f, Hx = Ht)
		nlp_stop.meta.nb_of_stop = iter
		calls = [nlp.counters.neval_obj,  nlp.counters.neval_grad, nlp.counters.neval_hess, nlp.counters.neval_hprod]
    end # while !OK


    xopt = xt
    fopt = ft
    ∇fopt = ∇f
    norm_∇fopt = norm_∇f

    # niter = verysucc + succ + unsucc;
    # verbose && @printf(" %10i unsuccessful iterations \n", unsucc)
    # verbose && @printf(" %10i successful iterations \n", succ)
    # verbose && @printf(" %10i very successful iterations \n", verysucc)
    # verbose && @printf(" %10i total successful iterations \n", verysucc + succ)
    # verbose && @printf(" %10i total iterations \n", niter)

    # calls = [nlp.counters.neval_obj, nlp.counters.neval_grad, nlp.counters.neval_hess, nlp.counters.neval_hprod]

    total_calls = nlp.counters.neval_obj + nlp.counters.neval_grad + n * nlp.counters.neval_hess + nlp.counters.neval_hprod
     # verbose && @printf(" %10i total f %10i g %10i H calls  %10i Hv calls  \n", nlp.counters.neval_obj, nlp.counters.neval_grad, nlp.counters.neval_hess, nlp.counters.neval_hprod)


    # if OK && verbose  @printf("\n     Final cost(x) = %f  ||g|| = %9.2e  nb of f/g/H evals = %i  converged\n", fopt, norm_∇fopt, total_calls)
    # elseif verbose
        # @printf("\n     Final cost(x) = %f  ||g|| = %9.2e  nb of f/g/H evals  = %i  NOT converged\n\n", fopt, norm_∇fopt, total_calls)
    # end;
    return nlp_at_x, nlp_stop.meta.optimal
end
