export TRARC

struct TRARCWorkspace
	xt
	xtnext
	d
	∇f
    ∇fnext
	function TRARCWorkspace(::Type{T}, ::Type{S}, n) where {T, S}
		return new(
			S(undef, n), # xt
			S(undef, n), # xtnext
			S(undef, n), # d
			S(undef, n), # ∇f
			S(undef, n), # ∇fnext
		)
	end
end

function TRARC(
	nlp 	:: AbstractNLPModel{T, S},
    nlp_stop	:: NLPStopping;
    TR 	:: TrustRegion = TrustRegion(T(10.0)),
	c 	:: Combi = Combi{T}(hessian_dense, PDataLDLt{T}, solve_modelTRDiag, preprocessLDLt, decreaseFact, Tparam{T}()),
	robust 	:: Bool = true,
    verbose  :: Bool = false,
) where {T, S}
	nlp_at_x = nlp_stop.current_state
    hessian_rep, PData, solve_model, pre_process, decrease, params = extract(c)
	workspace = TRARCWorkspace(T, S, nlp.meta.nvar)
	xt, xtnext, d, ∇f, ∇fnext = workspace.xt, workspace.xtnext, workspace.d, workspace.∇f, workspace.∇fnext

    α = TR.α₀  # initial Trust Region size
	max_unsuccinarow = TR.max_unsuccinarow

	xt .= nlp_at_x.x
    ft, ∇f = objgrad!(nlp, xt, ∇f)
	ftnext = ft
    norm_∇f = norm(∇f)
	nlp_stop.meta.optimality0 = norm_∇f

	OK = update_and_start!(nlp_stop, x = xt, fx = ft, gx = ∇f)
	!OK && Stopping.update!(nlp_at_x, Hx = hessian_rep(nlp, xt), convert = true)

	iter = 0 # counter different than stop count
    succ, unsucc, verysucc, unsuccinarow = 0, 0, 0, 0 # There is an issue with the way we cound unsuccinarow

	verbose && @info log_header([:iter, :f, :nrm_g, :λ,  :status, :α, :nrm_dtr, :nrm_dHO, :f_p_dTR, :f_p_dHO, :dir, :ΔqN], [Int64, T, T, T, String, T, T, T, T, T, String, T])
	verbose && @info log_row(Any[iter, nlp_at_x.fx, norm_∇f, 0.0, "First iteration", α])

	global dHO = nothing
	global dir = nothing # What is this information ?

    while !OK
		calls = [nlp.counters.neval_obj,  nlp.counters.neval_grad, nlp.counters.neval_hess, nlp.counters.neval_hprod]
        PData = pre_process(nlp_at_x.Hx, ∇f, params, calls, nlp_stop.meta.max_eval)

        if ~PData.OK
			@warn("Something wrong with PData")
			return nlp_at_x, nlp_stop.meta.optimal
		end

        success = false
        while !success & (unsuccinarow < max_unsuccinarow)
			d, dHO, λ = solve_model(nlp_stop, PData, α)
			dir = d == dHO ? "dTR == dHO" : "dTR != dHO" # Some unknown information?

            Δq = -(∇f + 0.5 * nlp_at_x.Hx * d)⋅d

            if Δq < 0.0 println("*******   Ascent direction in SolveModel: Δq = $Δq")
                println("  g⋅d = $(∇f⋅d), 0.5 d'Hd = $(0.5*(nlp_at_x.Hx*d)⋅d)  α = $α  λ = $λ")
                # cond issue with H?
				return nlp_at_x, nlp_stop.meta.optimal
            end
            slope = ∇f ⋅ d
			xt_original = copy(xt) # only used in prints, related to d et dHO
			xtnext .= xt .+ d

            iter += 1

			######################################
			# Why not just ...?
			# ftnext = obj(nlp, xtnext)
            try
                ftnext = obj(nlp, xtnext)
            catch
                ftnext = Inf;
            end

            if isnan(ftnext)
                ftnext = Inf
            end
			######################################
            Δf = ft - ftnext

            r, good_grad, ∇fnext = compute_r(nlp, ft, Δf, Δq, slope, d, xtnext, ∇fnext, robust)

            if r < TR.acceptance_threshold # unsucessful
				verbose && @info log_row(Any[iter, nlp_at_x.fx, norm_∇f, λ, "U", α, norm(d), norm(dHO), obj(nlp, xt_original .+ d), obj(nlp, xt_original .+ dHO), dir, Δq])
	        	unsucc += 1
	        	unsuccinarow += 1
	        	α = decrease(PData, α, TR)
	    	else
	        	success = true
	        	unsuccinarow = 0

	        	xt .= xtnext
	        	ft = ftnext
	        	if good_grad
					∇f .= ∇fnext
            	else
					grad!(nlp, xt, ∇f)
                end

				verysucc += 1
		        if r > TR.increase_threshold # very sucessful
		            α = increase(PData, α, TR)
					verbose && @info log_row(Any[iter, ft, stop_norm(∇f), λ, "V", α, norm(d), norm(dHO), obj(nlp, xt_original + d), obj(nlp, xt_original + dHO), dir, Δq])
		        else # sucessful
	                if r < TR.reduce_threshold
	                    α = decrease(PData, α, TR)
	                end
					verbose && @info log_row(Any[iter, ft, stop_norm(∇f), λ, "S", α, norm(d), norm(dHO),  obj(nlp, xt_original + d), obj(nlp, xt_original + dHO), dir, Δq])
		            succ += 1
		        end
	    	end
        end # while !success

		OK = update_and_stop!(nlp_stop, x = xt, fx = ft, gx = ∇f)
		success && Stopping.update!(nlp_at_x, Hx = hessian_rep(nlp, xt))
		nlp_stop.meta.nb_of_stop = iter
    end # while !OK

    return nlp_at_x, nlp_stop.meta.optimal
end
