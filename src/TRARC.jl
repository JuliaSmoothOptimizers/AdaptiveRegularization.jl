function TRARC(nlp :: AbstractNLPModel,
               x₀ :: Array{Float64,1},
               TR :: TrustRegion,
			   c :: Combi;
               atol :: Float64 = 1e-8,
               rtol :: Float64 = 1.0e-6,
               max_eval :: Int = 5000,
               max_iter :: Int=20000,
               #max_f :: Int=5000,
               #max_eval :: Int = 5000,
               robust :: Bool = true,
               verbose :: Bool = true
               )
	@show x₀
    hessian_rep,PData,solve_model,pre_process,decrease,params = extract(c)
	@show hessian_rep
	@show PData
	@show solve_model
	@show pre_process
	@show decrease
	@show params

    α = TR.α₀
    x, xnext, d, Df = copy(x₀), copy(x₀), copy(x₀), 0.0
    xopt = x
    λ = 1.0

    n = length(x)
	g = Vector{Float64}(undef, n)
	gnext = Vector{Float64}(undef, n)

    f = obj(nlp,x)
    if isnan(f) | (f==Inf)
        OK = false
	println("f nan or ∞");
	xopt=f
	gopt=f
	fopt=f
        niter = Inf
        calls = [Inf, Inf, Inf, Inf]
	return xopt,gopt,gopt,fopt,niter,calls,OK,:Bad_x0
    end

    fopt = f
    grad!(nlp,x,g)
    norm_g = stop_norm(g)
    norm_g0 = norm(g)
    gopt = g
    H = hessian_rep(nlp,x)

    optimal = (norm_g < atol) || (isinf(f) & (f<0.0))
    tired = false
    stalled = false
    finish = optimal || tired || stalled

    fnext = f
    iter = 0
    verbose && display_header_iterations()
    verbose && display_success(iter,fnext,norm_g0,0.0,α)

    succ, unsucc, verysucc, unsuccinarow = 0, 0, 0, 0

    max_unsuccinarow = TR.max_unsuccinarow

    calls = [0, 0, 0, 0]

    while ~finish

        PData = pre_process(H,g,params,calls,max_eval)

        if ~PData.OK return x,fopt,norm_g,Inf,Inf,[Inf,Inf,Inf,Inf],false,:PreFailed;  end

        success = false

#        while ~success & (iter < max_iter) & (unsuccinarow < TR.max_unsuccinarow)
        while ~success & (iter < max_iter) & (~stalled)
            try
                d,λ = solve_model(PData,α)
            catch
                println(" Problem in solve_model")

                return x,fopt,norm_g,norm(gopt),iter,calls,false,:AscentDir
            end

            Δq = -(g + 0.5 * H * d)⋅d

            if Δq < 0.0 println("*******   Ascent direction in SolveModel: Δq = $Δq")
                println("  g⋅d = $(g⋅d), 0.5 d'Hd = $(0.5*(H*d)⋅d)  α = $α  λ = $λ")
                println("  calls  $calls   iters $verysucc; $succ; $unsucc")

                return x,fopt,norm_g,norm(gopt),iter,calls,false,:AscentDir
            end
            slope = g ⋅ d
            xnext = x + d

            iter += 1

            # Trap illegal values for f, descent will backtrack until a legal value is reached
            try
                fnext = obj(nlp,xnext)
            catch
                fnext = Inf;
            end
            if isnan(fnext)
                fnext=Inf
            end

            Δf = f - fnext

            r, good_grad, gnext = compute_r(nlp,f,Δf,Δq,slope,d,xnext,gnext,robust)

            if r<TR.acceptance_threshold
                verbose && display_failure(iter,fnext,λ,α)
	        unsucc=unsucc+1
	        unsuccinarow = unsuccinarow +1
                stalled = unsuccinarow >= max_unsuccinarow
	        try
                    α = decrease(PData, α, TR)
                catch
                    stalled = true
                end
                #fbidon = obj(nlp,x)
	    else
	        success = true

	        unsuccinarow = 0
	        x = copy(xnext)
	        f = fnext

	        if good_grad  g = copy(gnext)
                else grad!(nlp,x,g);
                end

                norm_g=stop_norm(g)
                H = hessian_rep(nlp,x)
	        if r > TR.increase_threshold
	            α = increase(PData, α, TR)
	            verbose && display_v_success(iter,f,norm_g,λ,α)
                    verysucc += 1
	        else
                    if r < TR.reduce_threshold
                        α = decrease(PData, α, TR)
                    end
	            verbose && display_success(iter,f,norm_g,λ,α)
	            succ += 1
	        end
	    end
        end
        if isnan(f)  OK=false
	    println("f nan");
	    xopt=f
	    gopt=f
	    fopt=f
            niter = Inf
            calls = [Inf, Inf, Inf, Inf]
	    return xopt,gopt,fopt,niter,calls,OK
        end
        calls = [nlp.counters.neval_obj,  nlp.counters.neval_grad, nlp.counters.neval_hess, nlp.counters.neval_hprod]

        optimal = (norm_g < atol)| (norm_g <( rtol * norm_g0)) | (isinf(f) & (f<0.0))
        tired = (iter >= max_iter) | (sum(calls) > max_eval)



        finish = (optimal || tired || stalled)
    end

    xopt = x
    fopt = f
    gopt = g

    niter = verysucc+succ+unsucc;
    verbose && @printf(" %10i unsuccessful iterations \n",unsucc)
    verbose && @printf(" %10i successful iterations \n",succ)
    verbose && @printf(" %10i very successful iterations \n",verysucc)
    verbose && @printf(" %10i total successful iterations \n",verysucc+succ)
    verbose && @printf(" %10i total iterations \n",niter)

    OK = optimal

    calls = [nlp.counters.neval_obj, nlp.counters.neval_grad, nlp.counters.neval_hess, nlp.counters.neval_hprod]

    total_calls = nlp.counters.neval_obj + nlp.counters.neval_grad + n*nlp.counters.neval_hess+ nlp.counters.neval_hprod
     verbose && @printf(" %10i total f %10i g %10i H calls  %10i Hv calls  \n",nlp.counters.neval_obj, nlp.counters.neval_grad, nlp.counters.neval_hess, nlp.counters.neval_hprod)


    if OK && verbose
		@printf("\n     Final cost(x) = %f  ||g|| = %9.2e  number of f and g and H evals = %i  converged\n",fopt, norm(gopt),total_calls)
    elseif verbose
        @printf("\n     Final cost(x) = %f  ||g|| = %9.2e  number of f and g and H evals  = %i  NOT converged\n\n",fopt, norm(gopt),total_calls)
    end;

    if OK status = :Optimal
    elseif stalled status = :stalled
    else status = :UserLimit
    end

    return xopt,fopt,norm_g,norm(gopt),niter,calls,OK,status

end
