function preprocessST(H, g, params::TparamsST, n1, n2)
    τ = params.τ

    return  PDataST(H, g, τ, true)
end
