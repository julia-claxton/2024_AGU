import Pkg
include("$(@__DIR__)/code/Julia_ELFIN_Tools/install_dependencies.jl")
include("$(@__DIR__)/code/EPPBackscatterSimulation/install_dependencies.jl")
Pkg.add("BenchmarkTools")
Pkg.add("Profile")
Pkg.add("TickTock")