fname = :ARCSpectral

c = ARCTR.Combi(ARCTR.hessian_dense, ARCTR.PDataSpectral{Nothing}, ARCTR.solve_modelARCDiag, ARCTR.preprocessSpectral, ARCTR.decreaseFact, ARCTR.Tparam{Nothing}())

include("Template.jl")
