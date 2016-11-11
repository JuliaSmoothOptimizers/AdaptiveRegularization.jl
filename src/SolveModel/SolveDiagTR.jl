function solve_diagTR(λ,Δ,g̃,δ,ϵ)
    # λ underestimates the required value λstar such that ||d̃(λstar)|| = δ
    # where d̃(λ) =  -(Δ+λ*M) .\ g̃
    M = ones(Δ)
    λin = λ
    d̃ = -(Δ+λ*M) .\ g̃
    d̃d̃ = d̃⋅d̃
    normd̃ = sqrt(d̃d̃)
    
    tolerance = 1e-06
    # Newton iterations
    iter_nwt = 0
    #println(" Nwt $iter_nwt : λ = $λ  normd̃ = $normd̃  δ = $δ)")
    while (normd̃ >= (δ+tolerance*normd̃)) && (iter_nwt<40) 
        dotd̃ = (Δ+λ*M) .\ d̃
        Δλ = ((normd̃-δ)/δ) * (d̃d̃/(d̃ ⋅ dotd̃))
        λ =  max(λ + Δλ, λin+1.0e-10)
        d̃ = -(Δ+λ*M) .\ g̃
        d̃d̃ = d̃⋅d̃
        normd̃ = sqrt(d̃d̃)
        assert(λ>=0.0)
        iter_nwt += 1
    end
#println(" λ computation, iter_nwt : $iter_nwt  λ  $λ")

    return d̃,λ
end
