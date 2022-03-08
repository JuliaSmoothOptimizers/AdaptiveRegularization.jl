using ForwardDiff

import NLPModels.obj, NLPModels.grad, NLPModels.grad!, NLPModels.cons, NLPModels.cons!
import NLPModels.jac, NLPModels.jprod, NLPModels.jprod!, NLPModels.jtprod
import NLPModels.jtprod!, NLPModels.hess, NLPModels.hprod, NLPModels.hprod!
import NLPModels.increment!

export HOADNLPModel, obj, grad, grad!, cons, cons!, jac, jprod,
       jprod!, jtprod, jtprod!, hess, hprod, hprod!, ∇f³xuv, ∇f³xuv2

"""HOADNLPModel is an AbstractNLPModel using ForwardDiff to compute the
derivatives. HOHOADNLPModel uses ForwardDiff to compute 3rd order information.
This interface is very similar to HOADNLPModel with few additional functions which
are documented here.
For further documentation, refer to the HOADNLPModel.
"""
mutable struct HOADNLPModel <: AbstractNLPModel
  meta     :: NLPModelMeta

  counters :: Counters

  # Functions
  f :: Function
  c :: Function

  # AD Packages
  AD :: Any
end

function HOADNLPModel(f::Function, x0::AbstractVector; y0::AbstractVector = eltype(x0)[],
                      lvar::AbstractVector = eltype(x0)[], uvar::AbstractVector = eltype(x0)[],
                      lcon::AbstractVector = eltype(x0)[], ucon::AbstractVector = eltype(x0)[],
                      c::Function = (args...)->throw(NotImplementedError("cons")),
                      name::String = "Generic", lin::AbstractVector{Int}=Int[],
                      AD :: Any = ForwardDiff)

  nvar = length(x0)
  length(lvar) == 0 && (lvar = -Inf*ones(nvar))
  length(uvar) == 0 && (uvar =  Inf*ones(nvar))
  ncon = maximum([length(lcon); length(ucon); length(y0)])

  A = AD.hessian(f, x0)
  for i = 1:ncon
    A += AD.hessian(x->c(x)[i], x0) * (-1)^i
  end
  nnzh = nvar * (nvar + 1) / 2
  nnzj = 0

  if ncon > 0
    length(lcon) == 0 && (lcon = -Inf*ones(ncon))
    length(ucon) == 0 && (ucon =  Inf*ones(ncon))
    length(y0) == 0   && (y0 = zeros(ncon))
    nnzj = nvar * ncon
  end
  nln = setdiff(1:ncon, lin)

  meta = NLPModelMeta(nvar, x0=x0, lvar=lvar, uvar=uvar, ncon=ncon, y0=y0,
    lcon=lcon, ucon=ucon, nnzj=nnzj, nnzh=nnzh, lin=lin, nln=nln, minimize=true,
    islp=false, name=name)

  return HOADNLPModel(meta, Counters(), f, c, AD)
end

function obj(nlp :: HOADNLPModel, x :: AbstractVector)
  # increment!(nlp, :neval_obj)
  return nlp.f(x)
end

function grad(nlp :: HOADNLPModel, x :: AbstractVector)
  # increment!(nlp, :neval_grad)
  return ForwardDiff.gradient(nlp.f, x)
end

function grad!(nlp :: HOADNLPModel, x :: AbstractVector, g :: AbstractVector)
  # increment!(nlp, :neval_grad)
  ForwardDiff.gradient!(view(g, 1:length(x)), nlp.f, x)
  return g
end

function cons(nlp :: HOADNLPModel, x :: AbstractVector)
  # increment!(nlp, :neval_cons)
  return nlp.c(x)
end

function cons!(nlp :: HOADNLPModel, x :: AbstractVector, c :: AbstractVector)
  increment!(nlp, :neval_cons)
  c[1:nlp.meta.ncon] = nlp.c(x)
  return c
end

function jac(nlp :: HOADNLPModel, x :: AbstractVector)
  # increment!(nlp, :neval_jac)
  return ForwardDiff.jacobian(nlp.c, x)
end

function jac_coord(nlp :: HOADNLPModel, x :: AbstractVector)
  Jx = jac(nlp, x)
  m, n = nlp.meta.ncon, nlp.meta.nvar
  if VERSION < v"1.0"
    I = [(i,j) for i = 1:m, j = 1:n]
  else
    I = ((i,j) for i = 1:m, j = 1:n)
  end
  return (getindex.(I, 1)[:], getindex.(I, 2)[:], Jx[:])
end

function jprod(nlp :: HOADNLPModel, x :: AbstractVector, v :: AbstractVector)
  # increment!(nlp, :neval_jprod)
  return ForwardDiff.jacobian(nlp.c, x) * v
end

function jprod!(nlp :: HOADNLPModel, x :: AbstractVector, v :: AbstractVector, Jv :: AbstractVector)
  # increment!(nlp, :neval_jprod)
  Jv[1:nlp.meta.ncon] = ForwardDiff.jacobian(nlp.c, x) * v
  return Jv
end

function jtprod(nlp :: HOADNLPModel, x :: AbstractVector, v :: AbstractVector)
  # increment!(nlp, :neval_jtprod)
  return ForwardDiff.jacobian(nlp.c, x)' * v
end

function jtprod!(nlp :: HOADNLPModel, x :: AbstractVector, v :: AbstractVector, Jtv :: AbstractVector)
  # increment!(nlp, :neval_jtprod)
  Jtv[1:nlp.meta.nvar] = ForwardDiff.jacobian(nlp.c, x)' * v
  return Jtv
end

function hess(nlp :: HOADNLPModel, x :: AbstractVector; obj_weight :: Real = one(eltype(x)), y :: AbstractVector = eltype(x)[])
  increment!(nlp, :neval_hess)
  Hx = obj_weight == 0.0 ? spzeros(nlp.meta.nvar, nlp.meta.nvar) :
       ForwardDiff.hessian(nlp.f, x) * obj_weight
  for i = 1:min(length(y), nlp.meta.ncon)
    if y[i] != 0.0
      Hx += ForwardDiff.hessian(x->nlp.c(x)[i], x) * y[i]
    end
  end
  return tril(Hx)
end

function hess_coord(nlp :: HOADNLPModel, x :: AbstractVector; obj_weight :: Real = one(eltype(x)), y :: AbstractVector = eltype(x)[])
  Hx = hess(nlp, x, obj_weight=obj_weight, y=y)
  n = nlp.meta.nvar
  if VERSION < v"1.0"
    I = [(i,j,Hx[i,j]) for i = 1:n, j = 1:n if i ≥ j]
  else
    I = ((i,j,Hx[i,j]) for i = 1:n, j = 1:n if i ≥ j)
  end
  return (getindex.(I, 1), getindex.(I, 2), getindex.(I, 3))
end

function hprod(nlp :: HOADNLPModel, x :: AbstractVector, v :: AbstractVector;
               obj_weight :: Real = one(eltype(x)), y :: AbstractVector = eltype(x)[])
  Hv = zeros(nlp.meta.nvar)
  return hprod!(nlp, x, v, Hv, obj_weight=obj_weight, y=y)
end

function hprod!(nlp :: HOADNLPModel, x :: AbstractVector, v :: AbstractVector, Hv :: AbstractVector;
                obj_weight :: Real = one(eltype(x)), y :: AbstractVector = eltype(x)[])
  increment!(nlp, :neval_hprod)
  n = nlp.meta.nvar
  Hv[1:n] = obj_weight == 0.0 ? zeros(nlp.meta.nvar) :
          ForwardDiff.hessian(nlp.f, x) * v * obj_weight
  for i = 1:min(length(y), nlp.meta.ncon)
    if y[i] != 0.0
      Hv[1:n] += ForwardDiff.hessian(x->nlp.c(x)[i], x) * v * y[i]
    end
  end
  return Hv
end


function ∇f³xuv(nlp :: HOADNLPModel, x :: AbstractVector, u :: AbstractVector, v :: AbstractVector)
    ∇f(x) = ForwardDiff.gradient(nlp.f, x)
    # println("On  a ∇f")
    ϕᵤ(x) = ∇f(x)'*u
    # println("On  a ϕᵤ")
    ∇ϕᵤ(x) = ForwardDiff.gradient(ϕᵤ, x)
    # println("On  a ∇ϕᵤ")
    ψᵥ(x) =  ∇ϕᵤ(x)' * v
    # println("On  a ψᵥ")
    ∇ψᵥ(x) = ForwardDiff.gradient(ψᵥ, x)
    # println("On  a ∇ψᵥ")
    # increment!(nlp, :neval_Tuvprod)
    return ∇ψᵥ(x)
end

function ∇f³xuv2(nlp :: HOADNLPModel, x :: AbstractVector, u :: AbstractVector, v :: AbstractVector)
    ∇²fu(x) = ForwardDiff.hessian(nlp.f, x) * u
    ∇²fuv(x) = ∇²fu(x)' * v
    return ForwardDiff.gradient(∇²fuv, x)
end
