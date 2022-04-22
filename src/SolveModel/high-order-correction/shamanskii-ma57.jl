export Shamanskii_MA57

function Shamanskii_MA57(nlp_stop, PData::PDataFact, dₙ::Vector, gt::Vector)
    nlp_at_x = nlp_stop.current_state
    x = nlp_at_x.x
    xdemi = x + dₙ
    gtemp = grad(nlp_stop.pb, xdemi)
    ϵ2 = sqrt(eps(eltype(gt)))
    Γ = max.(abs.(PData.Δ), ϵ2)

    sg = PData.s .* gtemp
    ĝ = PData.L \ (sg[PData.pp])
    g̃ = PData.Q' * ĝ
    d̃ = -(PData.Δ) .\ g̃
    d = AInv(PData, d̃)

    dₕₒ = dₙ .+ d
    return dₕₒ
end
