if !isdir("../data") mkdir("../data") end
if !isdir("../data/horizons") mkdir("../data/horizons") end
println("build.jl running from ", pwd())

#include("build_julia_specific_conda.jl")
include("pyimports.jl")
