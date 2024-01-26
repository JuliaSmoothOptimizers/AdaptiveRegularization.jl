const solvers_const = Dict(
  :ARCqKdense =>
    (HessDense, PDataKARC, [(:shifts => 10.0 .^ (collect(-20.0:1.0:20.0)))]),
  :ARCqKOp =>
    (HessOp, PDataKARC, [:shifts => 10.0 .^ (collect(-20.0:1.0:20.0))]),
  :ARCqKsparse =>
    (HessSparse, PDataKARC, [:shifts => 10.0 .^ (collect(-20.0:1.0:20.0))]),
  :ARCqKCOO =>
    (HessSparseCOO, PDataKARC, [:shifts => 10.0 .^ (collect(-20.0:1.0:20.0))]),
  :ST_TRdense => (HessDense, PDataST, ()),
  :ST_TROp => (HessOp, PDataST, ()),
  :ST_TRsparse => (HessSparse, PDataST, ()),
  :TRKdense => (HessDense, PDataTRK, ()),
  :TRKOp => (HessOp, PDataTRK, ()),
  :TRKsparse => (HessSparse, PDataTRK, ()),
)

const solvers_nls_const = Dict(
  :ARCqKOpGN => (
    HessGaussNewtonOp,
    PDataKARC,
    [:shifts => 10.0 .^ (collect(-10.0:0.5:20.0))],
  ),
  :ST_TROpGN => (HessGaussNewtonOp, PDataST, ()),
  :ST_TROpGNLSCgls =>
    (HessGaussNewtonOp, PDataNLSST, [:solver_method => :cgls]),
  :ST_TROpGNLSLsqr =>
    (HessGaussNewtonOp, PDataNLSST, [:solver_method => :lsqr]),
  :ST_TROpLS => (HessOp, PDataNLSST, ()),
)
