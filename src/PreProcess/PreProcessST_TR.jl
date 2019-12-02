function preprocessST(H, g, params::TparamsST, n1, n2)
    ζ = params.ζ

    return  PDataST(H, g, ζ, true)
end
