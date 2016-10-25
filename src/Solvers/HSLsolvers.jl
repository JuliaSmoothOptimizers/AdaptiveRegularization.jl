using HSL

HSL_solvers = Function[]


include("ARCMA97.jl")
push!(HSL_solvers,eval(fname))

include("ARCMA97_abs.jl")
push!(HSL_solvers,eval(fname))

include("TRMA97.jl")
push!(HSL_solvers,eval(fname))

include("TRMA97_abs.jl")
push!(HSL_solvers,eval(fname))

