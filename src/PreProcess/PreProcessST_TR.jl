function preprocessST(PData::PDataST, H, g, params::TparamsST, n1, n2)
    ζ = params.ζ
    PData.ζ = ζ
    PData.g = g
    PData.H = H
    PData.OK = true

    return PData
end
