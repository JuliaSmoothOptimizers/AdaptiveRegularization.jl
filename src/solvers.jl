solvers_const = Dict(
    :ARCLDLt_abs => (HessDense, PDataLDLt, solve_modelARCDiagAbs, ()),
    :ARCLDLt => (HessDense, PDataLDLt, solve_modelARCDiag, ()),
    # :ARCMA57_abs => (HessSparse, PDataMA57, solve_modelARCDiagAbs, ()),
    # :ARCMA57_Sham_vs_Nwt => (HessSparse, PDataMA57, solve_modelARCDiag, ()),
    # :ARCMA57 => (HessSparse, PDataMA57, solve_modelARCDiag, ()),
    # :ARCMA97_abs => (HessDense, PDataMA97, solve_modelARCDiagAbs, ()),
    # :ARCMA97 => (HessDense, PDataMA97, solve_modelARCDiag, ()),
    :ARCqKdense => (
        HessDense,
        PDataKARC,
        solve_modelKARC,
        [(:shifts => 10.0 .^ (collect(-20.0:1.0:20.0)))],
    ),
    :ARCqKOp => (
        HessOp,
        PDataKARC,
        solve_modelKARC,
        [:shifts => 10.0 .^ (collect(-20.0:1.0:20.0))],
    ),
    :ARCqKsparse => (
        HessSparse,
        PDataKARC,
        solve_modelKARC,
        [:shifts => 10.0 .^ (collect(-20.0:1.0:20.0))],
    ),
    :ARCSpectral_abs => (HessDense, PDataSpectral, solve_modelARCDiagAbs, ()),
    :ARCSpectral => (HessDense, PDataSpectral, solve_modelARCDiag, ()),
    :ST_TRdense => (HessDense, PDataST, solve_modelST_TR, ()),
    :ST_TROp => (HessOp, PDataST, solve_modelST_TR, ()),
    :ST_TRsparse => (HessSparse, PDataST, solve_modelST_TR, ()),
    :TRKdense => (HessDense, PDataKTR, solve_modelKTR, ()),
    :TRKOp => (HessOp, PDataKTR, solve_modelKTR, ()),
    :TRKsparse => (HessSparse, PDataKTR, solve_modelKTR, ()),
    :TRLDLt_abs => (HessDense, PDataLDLt, solve_modelTRDiagAbs, ()),
    :TRLDLt => (HessDense, PDataLDLt, solve_modelTRDiag, ()),
    # :TRMA57_abs => (HessSparse, PDataMA57, solve_modelTRDiagAbs, ()),
    # :TRMA57 => (HessSparse, PDataMA57, solve_modelTRDiag, ()),
    # :TRMA97_abs => (HessDense, PDataMA97, solve_modelTRDiagAbs, ()),
    # :TRMA97 => (HessSparse, PDataMA97, solve_modelTRDiag, ()),
    :TRSpectral_abs => (HessDense, PDataSpectral, solve_modelTRDiagAbs, ()),
    :TRSpectral => (HessDense, PDataSpectral, solve_modelTRDiag, ()),
)
