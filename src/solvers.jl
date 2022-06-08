const solvers_const = Dict(
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
    :ARCqKCOO => (
        HessSparseCOO,
        PDataKARC,
        solve_modelKARC,
        [:shifts => 10.0 .^ (collect(-20.0:1.0:20.0))],
    ),
    :ST_TRdense => (HessDense, PDataST, solve_modelST_TR, ()),
    :ST_TROp => (HessOp, PDataST, solve_modelST_TR, ()),
    :ST_TRsparse => (HessSparse, PDataST, solve_modelST_TR, ()),
    :TRKdense => (HessDense, PDataKTR, solve_modelKTR, ()),
    :TRKOp => (HessOp, PDataKTR, solve_modelKTR, ()),
    :TRKsparse => (HessSparse, PDataKTR, solve_modelKTR, ()),
)
