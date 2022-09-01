var documenterSearchIndex = {"docs":
[{"location":"reference/#Reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"​","category":"page"},{"location":"reference/#Contents","page":"Reference","title":"Contents","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"​","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"Pages = [\"reference.md\"]","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"​","category":"page"},{"location":"reference/#Index","page":"Reference","title":"Index","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"​","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"Pages = [\"reference.md\"]","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"​","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"Modules = [ARCTR]","category":"page"},{"location":"reference/#ARCTR.HessDense","page":"Reference","title":"ARCTR.HessDense","text":"HessDense(::AbstractNLPModel{T,S}, n)\n\nReturn a structure used for the evaluation of dense Hessian matrix.\n\n\n\n\n\n","category":"type"},{"location":"reference/#ARCTR.HessOp","page":"Reference","title":"ARCTR.HessOp","text":"HessOp(::AbstractNLPModel{T,S}, n)\n\nReturn a structure used for the evaluation of the Hessian matrix as an operator.\n\n\n\n\n\n","category":"type"},{"location":"reference/#ARCTR.HessSparse","page":"Reference","title":"ARCTR.HessSparse","text":"HessSparse(::AbstractNLPModel{T,S}, n)\n\nReturn a structure used for the evaluation of sparse Hessian matrix.\n\n\n\n\n\n","category":"type"},{"location":"reference/#ARCTR.HessSparseCOO","page":"Reference","title":"ARCTR.HessSparseCOO","text":"HessSparseCOO(::AbstractNLPModel{T,S}, n)\n\nReturn a structure used for the evaluation of sparse Hessian matrix in COO-format.\n\n\n\n\n\n","category":"type"},{"location":"reference/#ARCTR.PDataKARC","page":"Reference","title":"ARCTR.PDataKARC","text":"PDataKARC(::Type{S}, ::Type{T}, n)\n\nReturn a structure used for the preprocessing of ARCqK methods.\n\n\n\n\n\n","category":"type"},{"location":"reference/#ARCTR.PDataKTR","page":"Reference","title":"ARCTR.PDataKTR","text":"PDataKTR(::Type{S}, ::Type{T}, n)\n\nReturn a structure used for the preprocessing of TRK methods.\n\n\n\n\n\n","category":"type"},{"location":"reference/#ARCTR.PDataST","page":"Reference","title":"ARCTR.PDataST","text":"PDataST(::Type{S}, ::Type{T}, n)\n\nReturn a structure used for the preprocessing of Steihaug-Toint methods.\n\n\n\n\n\n","category":"type"},{"location":"reference/#ARCTR.TRARCWorkspace","page":"Reference","title":"ARCTR.TRARCWorkspace","text":"TRARCWorkspace(nlp, ::Type{Hess}, n)\n\nPre-allocate the memory used during the TRARC call for the problem nlp of size n. The possible values for Hess are: HessDense, HessSparse, HessSparseCOO, HessOp. Return a TRARCWorkspace structure.\n\n\n\n\n\n","category":"type"},{"location":"reference/#ARCTR.TrustRegion","page":"Reference","title":"ARCTR.TrustRegion","text":"TrustRegion(α₀::T;kwargs...)\n\nSelect the main parameters used in the TRARC algorithm with α₀ as initial TR/ARC parameter. The keyword arguments are:\n\nmax_α::T: Maximum value for α. Default T(1) / sqrt(eps(T)).\nacceptance_threshold::T: Ratio over which the step is successful. Default T(0.1).\nincrease_threshold::T: Ratio over which we increase α. Default T(0.75).\nreduce_threshold::T: Ratio under which we decrease α. Default T(0.1).\nincrease_factor::T: Factor of increase of α. Default T(5.0).\ndecrease_factor::T: Factor of decrease of α. Default T(0.1).\nmax_unsuccinarow::Int: Limit on the number of successive unsucessful iterations. Default 30.\n\nReturns a TrustRegion structure.\n\nThis can be compared to https://github.com/JuliaSmoothOptimizers/SolverTools.jl/blob/main/src/trust-region/basic-trust-region.jl\n\n\n\n\n\n","category":"type"},{"location":"reference/#ARCTR.TrustRegionException","page":"Reference","title":"ARCTR.TrustRegionException","text":"Exception type raised in case of error.\n\n\n\n\n\n","category":"type"},{"location":"reference/#ARCTR.TRARC","page":"Reference","title":"ARCTR.TRARC","text":"TRARC(nlp; kwargs...)\n\nCompute a local minimum of an unconstrained optimization problem using trust-region (TR)/adaptive regularization with cubics (ARC) methods.\n\nArguments\n\nnlp::AbstractNLPModel: the model solved, see NLPModels.jl.\n\nThe keyword arguments include\n\nTR::TrustRegion: structure with trust-region/ARC parameters, see TrustRegion. Default: TrustRegion(T(10.0)).\nhess_type::Type{Hess}: Structure used to handle the hessian. The possible values are: HessDense, HessSparse, HessSparseCOO, HessOp. Default: HessOp.\npdata_type::Type{ParamData} Structure used for the preprocessing step. Default: PDataKARC.\nsolve_model::Function Function used to solve the subproblem. Default: solve_modelKARC.\nrobust::Bool: true implements a robust evaluation of the model. Default: true.\nverbose::Bool: true prints iteration information. Default: false.\n\nAdditional kwargs are used for stopping criterion, see Stopping.jl.\n\nOutput\n\nThe returned value is a GenericExecutionStats, see SolverCore.jl.\n\nThis implementation uses Stopping.jl. Therefore, it is also possible to used\n\nTRARC(stp; kwargs...)\n\nwhich returns the stp::NLPStopping updated.\n\nFor advanced usage, the principal call to the solver uses a TRARCWorkspace.\n\nTRARC(stp, pdata, workspace, trust_region_parameters; kwargs...)\n\nSome variants of TRARC are already implemented and listed in ARCTR.ALL_solvers.\n\nReferences\n\nThis method unifies the implementation of trust-region and adaptive regularization with cubics as described in\n\nDussault, J.-P. (2020).\nA unified efficient implementation of trust-region type algorithms for unconstrained optimization.\nINFOR: Information Systems and Operational Research, 58(2), 290-309.\n10.1080/03155986.2019.1624490\n\nExamples\n\njulia> using ARCTR, ADNLPModels\njulia> nlp = ADNLPModel(x -> 100 * (x[2] - x[1]^2)^2 + (x[1] - 1)^2, [-1.2; 1.0]);\njulia> stats = TRARC(nlp)\n\"Execution stats: first-order stationary\"\n\n\n\n\n\n","category":"function"},{"location":"reference/#ARCTR.compute_r-Union{Tuple{T}, Tuple{Any, T, Any, Any, Any, Any, Any, Any, Any}} where T","page":"Reference","title":"ARCTR.compute_r","text":"compute_r(nlp, f, Δf, Δq, slope, d, xnext, gnext, robust)\n\nCompute the actual vs predicted reduction ratio ∆f/Δq.\n\nArguments:\n\nnlp: Current model we are trying to solve\nf: current objective value\nΔf: = f - f_trial is the actual reduction is an objective/merit/penalty function,\nΔq: q - q_trial is the reduction predicted by the model q of f.\nslope: current slope\nd: potential next direction\nxnext: potential next iterate\ngnext: current gradient value, if good_grad is true, then this value has been udpated.\nrobust: if true, try to trap potential cancellation errors\n\nOutput:\n\nr: reduction ratio ∆f/Δq\ngood_grad: true if gnext has been recomputed\ngnext: gradient.\n\nWe assume that qis being minimized, and therefore thatΔq > 0`.\n\n\n\n\n\n","category":"method"},{"location":"reference/#ARCTR.decrease-Union{Tuple{T}, Tuple{ARCTR.TPData, T, TrustRegion}} where T","page":"Reference","title":"ARCTR.decrease","text":"decrease(X::TPData, α::T, TR::TrustRegion)\n\nReturn a decreased α.\n\n\n\n\n\n","category":"method"},{"location":"reference/#ARCTR.hessian!","page":"Reference","title":"ARCTR.hessian!","text":"hessian!(workspace::HessDense, nlp, x)\n\nReturn the Hessian matrix of nlp at x in-place with memory update of workspace.\n\n\n\n\n\n","category":"function"},{"location":"reference/#ARCTR.increase-Union{Tuple{T}, Tuple{ARCTR.TPData, T, TrustRegion}} where T","page":"Reference","title":"ARCTR.increase","text":"increase(X::TPData, α::T, TR::TrustRegion)\n\nReturn an increased α.\n\n\n\n\n\n","category":"method"},{"location":"reference/#ARCTR.preprocess-Tuple{ARCTR.TPData, Any, Any, Any, Any, Any, Any}","page":"Reference","title":"ARCTR.preprocess","text":"preprocess(PData::TPData, H, g, gNorm2, n1, n2, α)\n\nFunction called in the TRARC algorithm every time a new iterate has been accepted.\n\nArguments\n\nPData::TPData: data structure used for preprocessing.\nH: current Hessian matrix.\ng: current gradient.\ngNorm2: 2-norm of the gradient.\nn1: Current count on the number of Hessian-vector products.\nn2: Maximum number of Hessian-vector products accepted.\nα: current value of the TR/ARC parameter.\n\nIt returns PData.\n\n\n\n\n\n","category":"method"},{"location":"reference/#ARCTR.solve_model-Tuple{Any, Any, Any, Any, ARCTR.TPData, Any}","page":"Reference","title":"ARCTR.solve_model","text":"solve_model(H, g, gNorm2, nlp_stop, PData::TPData, α)\n\nFunction called in the TRARC algorithm to solve the subproblem.\n\nArguments\n\nPData::TPData: data structure used for preprocessing.\nH: current Hessian matrix.\ng: current gradient.\ngNorm2: 2-norm of the gradient.\nnlp_stop: Current NLPStopping representing the problem.\nα: current value of the TR/ARC parameter.\n\nIt returns a couple (PData.d, PData.λ). Current implementations include: solve_modelKARC, solve_modelKTR, solve_modelST_TR.\n\n\n\n\n\n","category":"method"},{"location":"benchmark/#Benchmarks","page":"Tutorial","title":"Benchmarks","text":"","category":"section"},{"location":"benchmark/#CUTEst-benchmark","page":"Tutorial","title":"CUTEst benchmark","text":"","category":"section"},{"location":"#ARCTR","page":"Introduction","title":"ARCTR","text":"","category":"section"}]
}