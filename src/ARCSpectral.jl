fname = :ARCSpectral
c = JPDJPD.Combi(JPDJPD.hessian_dense,JPDJPD.PDataSpectral,JPDJPD.solve_modelARCDiag,JPDJPD.preprocessSpectral,JPDJPD.decreaseFact,JPDJPD.Tparam())
@eval begin
    export $fname


    include("Template.jl")
end
