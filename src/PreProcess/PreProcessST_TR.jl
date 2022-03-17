function preprocess(PData::PDataST, H, g, n1, n2)
    ζ = PData.ζ
    PData.ζ = ζ
    PData.g = g
    PData.H = H
    PData.OK = true

    return PData
end
