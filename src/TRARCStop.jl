#@debug 
function TRARC(nlp :: AbstractNLPModel,
               x₀ :: Array{Float64,1},
               TR :: TrustRegion, c :: Combi;
               s :: TStopping = TStopping(nlp),
               verbose :: Bool = true
               )                

    hessian_rep,PData,solve_model,pre_process,decrease,params = extract(c)

    StopNorm = s.optimality_residual
     
    s = start!(s,x₀)

    α = TR.α₀
    x, xnext, d, Df = copy(x₀), copy(x₀), copy(x₀), 0.0
    xopt = x
    λ = 1.0
    
    n = length(x)
    ∇f = Array(Float64,n)
    ∇fnext = Array(Float64,n)    
    
    f = obj(nlp,x)
    if isnan(f) | (f==Inf)
        OK = false
	println("f nan or ∞");
	xopt=f
	∇fopt=f
	fopt=f
        niter = Inf
        calls = [Inf, Inf, Inf, Inf]
	return xopt,∇fopt,∇fopt,fopt,niter,calls,OK,:Bad_x0
    end

    fopt = f
    grad!(nlp,x,∇f)
    norm_∇f = stop_norm(∇f)
    norm_∇f0 = norm_∇f
    ∇fopt = ∇f
    norm_∇fopt = norm_∇f
    H = hessian_rep(nlp,x)
    
    #optimal = (norm_g < atol) || (isinf(f) & (f<0.0))
    #tired = false
    
    fnext = f
    iter = 0

    optimal, unbounded, tired, elapsed_time = stop(s,iter,x,f,∇f)
    stalled = false
    finish = optimal || unbounded || tired || stalled

    verbose && display_header_iterations()
    verbose && display_success(iter,fnext,norm_g0,0.0,α)
    
    succ, unsucc, verysucc, unsuccinarow = 0, 0, 0, 0
    
    max_unsuccinarow = TR.max_unsuccinarow
    
    calls = [0, 0, 0, 0]
    
    while ~finish

        PData = pre_process(H,∇f,params,calls,s.max_eval)
#                println("===>>> cond(H) = $(cond(H))  cond(L) = $(cond(PData.L))")
#                println(" ||H - P'LDL'P|| = $(norm(H - PData.P'*PData.L*PData.D*PData.L'*PData.P))")
        #if cond(PData.L) > 10000.0*cond(H) println("===>>> cond(H) = $(cond(H))  cond(L) = $(cond(PData.L))")
        #end
        
        if ~PData.OK return x,fopt,norm_∇f,Inf,Inf,[Inf,Inf,Inf,Inf],false,:PreFailed;  end
        
        success = false
        
        while ~success & (iter < s.max_iter) & (unsuccinarow < TR.max_unsuccinarow)
            try
                d,λ = solve_model(PData,α)
            catch
                println(" Problem in solve_model")

                return x,fopt,norm_∇f,norm_∇fopt,iter,calls,false,:AscentDir
            end
            #println("******* TRARC:  g⋅d = $(g⋅d), 0.5 d'*H*d = $(0.5*(H * d)⋅d)")

            Δq = -(∇f + 0.5*H*d)⋅d

            if Δq<0.0 println("*******   Ascent direction in SolveModel: Δq = $Δq")
                println("  g⋅d = $(∇f⋅d), 0.5 d'Hd = $(0.5*(H*d)⋅d)  α = $α  λ = $λ")
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

                return x,fopt,norm_∇f,norm_∇fopt,iter,calls,false,:AscentDir
            end    
            slope = ∇f ⋅ d
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
                        
            r, good_grad, ∇fnext = compute_r(nlp,f,Δf,Δq,slope,d,xnext,∇fnext)
            
            if r<TR.acceptance_threshold
                verbose && display_failure(iter,fnext,λ,α)
	        unsucc=unsucc+1
	        unsuccinarow = unsuccinarow +1
	        α = decrease(PData, α, TR)
                fbidon = obj(nlp,x)
	    else
	        success = true
                
	        unsuccinarow = 0
	        x = copy(xnext)
	        f = fnext
                
	        if good_grad  ∇f = copy(∇fnext)
                else grad!(nlp,x,∇f);
                end
                
                norm_∇f=stop_norm(∇f)
                H = hessian_rep(nlp,x)
	        if r > TR.increase_threshold
	            α = increase(PData, α, TR)
	            verbose && display_v_success(iter,f,norm_∇f,λ,α)
                    verysucc += 1
	        else
                    if r < TR.reduce_threshold
                        α = decrease(PData, α, TR)
                    end
	            verbose && display_success(iter,f,norm_∇f,λ,α)
	            succ += 1
	        end
	    end
        end
        if isnan(f)  OK=false 
	    println("f nan");
	    xopt=f
	    ∇fopt=f
	    fopt=f
            niter = Inf
            calls = [Inf, Inf, Inf, Inf]
	    return xopt,∇fopt,∇fopt,fopt,niter,calls,OK
        end
        calls = [nlp.counters.neval_obj,  nlp.counters.neval_grad, nlp.counters.neval_hess, nlp.counters.neval_hprod]

        optimal, unbounded, tired, elapsed_time = stop(s,iter,x,f,∇f)
        #optimal = (norm_g < atol)| (norm_g <( rtol * norm_g0)) | (isinf(f) & (f<0.0))
        #tired = (iter >= itmax) | (sum(calls) > max_calls)
        stalled = unsuccinarow >= max_unsuccinarow


        finish = optimal || unbounded || tired || stalled
    end
    
    xopt = x
    fopt = f
    ∇fopt = ∇f
    norm_∇fopt = norm_∇f

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
    
    
    if OK && verbose  @printf("\n     Final cost(x) = %f  ||g|| = %9.2e  \#f/g/H evals = %i  converged\n",fopt, norm_∇fopt,total_calls)
    elseif verbose
        @printf("\n     Final cost(x) = %f  ||g|| = %9.2e  \#f/g/H evals  = %i  NOT converged\n\n",fopt, norm_∇fopt,total_calls)
    end;  
    
    if     optimal   status = :Optimal
    elseif stalled   status = :Stalled
    elseif unbounded status = :Unbounded
    else             status = :UserLimit
    end
    
    return xopt,fopt,norm_∇fopt,norm(∇fopt),niter,calls,OK,status
    
end

       
