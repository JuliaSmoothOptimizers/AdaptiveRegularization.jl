using Pkg;
Pkg.activate("");
using CUTEst

nmax = 1000000 # maximum is 1000000
_pnames = CUTEst.select(contype = "unc", max_var = nmax)

# Remove all the problems ending by N.
pnamesNE = _pnames[findall(x -> occursin(r"NE\b", x), _pnames)]
pnames = setdiff(_pnames, pnamesNE)

open("list_problems_$nmax.dat", "w") do io
  for name in pnames
    write(io, name * "\n")
  end
end
