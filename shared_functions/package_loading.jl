#This script is used to load all of the packages required for the dKMC simulations.

#List of required packages.
required_packages = ["Dates", "DelimitedFiles", "Distributed", "Distributions", "LinearAlgebra", "Logging", "NumericalIntegration", "Measurements", "Polynomials", "QuadGK", "Random", "Statistics", "StatsBase", "TimerOutputs", "YAML"]

#Check if the packages have already been added, and if not add them.
using Pkg
for package in required_packages
    if Base.find_package(package) === nothing
        Pkg.add(package)
    end
end

#Load the required packages.
using Dates
using DelimitedFiles
using Distributed
using Distributions
using LinearAlgebra
using Logging
using NumericalIntegration
using Measurements
using Polynomials
using QuadGK
using Random
using Statistics
using StatsBase 
using TimerOutputs
using YAML