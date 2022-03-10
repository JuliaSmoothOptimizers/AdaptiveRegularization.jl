using HSL

T = Float64
# definitions applicable to all packages
const data_map = Dict{Type, Type}(Float32 => Cfloat,
                                  Float64 => Cdouble,
                                  Complex64 => Cfloat,
                                  Complex128 => Cdouble)

rows = Int32[1, 1, 2, 2, 3, 3, 5]
cols = Int32[1, 2, 3, 5, 3, 4, 5]
vals = T[2, 3, 4, 6, 1, 5, 1]
A = sparse(rows, cols, vals) ; A = A + triu(A, 1)'

H57 = convert(SparseMatrixCSC{data_map[T],Int32}, A)
L = tril(H57)
M = ma57_coord(L.n, findnz(L)...)


b = T[8, 45, 31, 15, 17]
xexact = T[1, 2, 3, 4, 5]

ϵ = sqrt(eps(eltype(A)))
ma57_factorize(M)
x = ma57_solve(M, b)

(L, D57, s, p) = ma57_get_factors(M)
S = spdiagm(s)
P = speye(M.n)[p, :]

r0 = P * S * A * S * P' - L * D57 * L'
pi = invperm(p)
r0b = S[p,:] * A * S[:,pi] - L * D57 * L'

Sm1 = spdiagm(1.0 ./ s)

x1 = ma57_solve(M, b, job = :LS)
r1 = (L*S * x1) - (P*S*b)
r1b = (L * x1) - (P*S*b)
x2 = ma57_solve(M, b, job = :LPS)
r2 = (L'*P*Sm1 * x2) - (Sm1*b)
r2b = (L'*P*Sm1 * x2) - b


H97 = convert(SparseMatrixCSC{data_map[T],Int64}, A)
PLD = Ma97(H97, print_level=-1)
ma97_factorize(PLD)

ret = ma97_inquire(PLD)
vD = ret[2]
pivot = ret[1]
# Take care of zero valued pivots
z(a)= a==0
pivot[find(z,pivot)] = collect(1:length(pivot))[find(z,pivot)]

vD1 = copy(vD[1:1,:])'[:,1:1]
vD2 = copy(vD[2:2,1:end-1])'[:,1:1]

D97 = spdiagm((vD1,vD2,vD2),[0 -1 1]) # don't forget inquire returns D^-1
