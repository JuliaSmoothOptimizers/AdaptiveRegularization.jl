#export ARCSpectral, ARCSpectral_abs, ARCLDLt, ARCqKOp
export ALL_solvers

# Valid combinations
#

ALL_solvers = Function[]


include("ARCSpectral.jl")      #works in 0.7/1.1
push!(ALL_solvers,eval(fname))

include("ARCSpectral_abs.jl") #works in 0.7/1.1
push!(ALL_solvers,eval(fname))



include("TRSpectral.jl") #works in 0.7/1.1
push!(ALL_solvers,eval(fname))

include("TRSpectral_abs.jl") #works in 0.7/1.1
push!(ALL_solvers,eval(fname))



include("ARCLDLt.jl")         #works in 0.7/1.1
push!(ALL_solvers,eval(fname))

include("ARCLDLt_abs.jl")      #works in 0.7/1.1
push!(ALL_solvers,eval(fname))



include("TRLDLt.jl")             #works in 0.7/1.1
push!(ALL_solvers,eval(fname))

include("TRLDLt_abs.jl")         #works in 0.7/1.1
push!(ALL_solvers,eval(fname))



include("ARCqKOp.jl")			   #works in 0.7/1.1
push!(ALL_solvers,eval(fname))


include("ARCqKsparse.jl")          #works in 0.7/1.1
push!(ALL_solvers,eval(fname))

include("ARCqKdense.jl")          #works in 0.7/1.1
push!(ALL_solvers,eval(fname))



include("TRKOp.jl")              #works in 0.7/1.1
push!(ALL_solvers,eval(fname))

include("TRKsparse.jl")          #works in 0.7/1.1
push!(ALL_solvers,eval(fname))

include("TRKdense.jl")           #works in 0.7/1.1
push!(ALL_solvers,eval(fname))



include("ST_TROp.jl")            #works in 0.7/1.1
push!(ALL_solvers,eval(fname))

include("ST_TRsparse.jl")        #works in 0.7/1.1
push!(ALL_solvers,eval(fname))

include("ST_TRdense.jl")        #works in 0.7/1.1
push!(ALL_solvers,eval(fname))



include("ST_ARCOp.jl")            #works in 0.7/1.1
push!(ALL_solvers,eval(fname))

include("ST_ARCsparse.jl")        #works in 0.7/1.1
push!(ALL_solvers,eval(fname))

include("ST_ARCdense.jl")         #works in 0.7/1.1
push!(ALL_solvers,eval(fname))


## Needs further work to understand the HSL methods
## Will update them in the future
# include("ARCMA97.jl")
# push!(ALL_solvers,eval(fname))
#
# include("ARCMA97_abs.jl")
# push!(ALL_solvers,eval(fname))
#
# include("TRMA97.jl")
# push!(ALL_solvers,eval(fname))
#
# include("TRMA97_abs.jl")
# push!(ALL_solvers,eval(fname))
#
#
#
#
# include("ARCMA57.jl")
# push!(ALL_solvers,eval(fname))
#
# include("ARCMA57_abs.jl")
# push!(ALL_solvers,eval(fname))
#
# include("TRMA57.jl")
# push!(ALL_solvers,eval(fname))
#
# include("TRMA57_abs.jl")
# push!(ALL_solvers,eval(fname))
