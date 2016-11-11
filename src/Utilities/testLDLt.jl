using Base.Test


n = 10;
A = rand(n, n); A = A + A';

include("ldlt_symm.jl")

L = Array(Float64,2)
D = Array(Float64,2)
P = Array(Float64,2)
pp = Array(Int,1)
ρ = Float64
ncomp = Int64

try (L, D, pp, ρ, ncomp) = ldlt_symm(A)
catch except
    @test rethrow(except)
end 

#@test_approx_eq P*A*P'  L*D*L'
@test_approx_eq A[pp,pp]  L*D*L'

v = rand(n,1);
#@test_approx_eq P*v  v[pp]
#@test_approx_eq P'*v  v[invperm(pp)]

