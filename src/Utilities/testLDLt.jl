using Base.Test


n = 10;
A = rand(n, n); A = A + A';

include("ldlt_symm.jl")

L = Array(Float64,2)
D = Array(Float64,2)
P = Array(Float64,2)
ρ = Float64
ncomp = Int64

try (L, D, P, ρ, ncomp) = ldlt_symm(A)
catch except
    @test rethrow(except)
end 

@test_approx_eq P*A*P'  L*D*L'
