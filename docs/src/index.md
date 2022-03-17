# ARCTR

We list here the possible usage options of the unified implementation.

## Preprocess

This step is called every time a new iteration is accepted, i.e. a successful or very successful step is accepted.
The package contains several implementations:
- preprocessKARC (PreprocessKARC.jl) return a `PDataK` and extra function: decreaseKARC(X::PDataK, α::Float64, TR::TrustRegion)
- preprocessKTR (PreprocessKTR.jl) return a `PDataK` and extra function: decreaseKTR(X::PDataK, α::Float64, TR::TrustRegion)
- preprocessLDLt (PreprocessLDLt.jl) return a `PDataLDLt` and extra function: AInv(X::PDataLDLt, d̃::Array{T,1}) where {T}, reconstructH(X::PDataLDLt)
- preprocessLDLt2 (PreprocessLDLt.jl) return a `PDataLDLt` and extra function: AInv(X::PDataLDLt, d̃::Array{T,1}) where {T}, reconstructH(X::PDataLDLt)
- preprocessMA57 (PreprocessMA57.jl) return a `PDataMA57` and extra function: AInv(X::PDataMA57, d̃::Array{T,1}) where {T}, reconstructH(X::PDataMA57)
- preprocessMA97 (PreprocessMA97.jl) return a `PDataMA97` and extra function: AInv(X::PDataMA97, d̃::Array{T,1}) where {T}, reconstructH(X::PDataMA97)
- preprocessSpectral (PreprocessSpectral.jl) return a `PDataSpectral` and extra function: AInv(X::PDataSpectral, d̃::Array{Float64,1}), reconstructH(X::PDataSpectral)
- preprocessST (PreprocessST_TR.jl) return a `PDataST` and extra function: none.

`AInv` is the function used to solve the linear system, I think.

The **decrease** functions ((X::TPData, α::T, TR::TrustRegion) where {T} -> T) are described in Utilies.jl and contains
- decreaseGen
- decreaseFact

**TODO**: 
-[ ] update for in-place preprocess functions
-[X] update for in-place `hessian_rep` functions
-[ ] we should probably move decreaseKTR and decreaseKARC here in Utilies.jl
-[ ] Why is there a LDLt and LDLt2? LDLt2 that depend on a `bunchkaufman` function (from LinearAlgebra?)
-[ ] we could use multiple-dispatch and have just one function `decrease` and `increase`.
-[X] remove `reconstructH` as it is nowhere used.
-[ ] extract parameters in KARC and KTR
-[ ] reuse KrylovSolver structure

## Parameter tuning

It seems some problems are not solved because of some parameter tuning:
- NELSONLS
