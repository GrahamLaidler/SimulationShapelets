using XLSX,
    Random,
    LaTeXStrings,
    #KernelDensity,
    Plots.Measures,
    JLD2,
    StatsPlots,
    LinearAlgebra,
    StatsBase,
    #InvertedIndices,
    Distributions,
    ProgressMeter

abstract type LocInvMethod end
struct Yes <: LocInvMethod end
struct No <: LocInvMethod end

abstract type SearchOver end
struct All <: SearchOver end
struct Class1 <: SearchOver end
struct Class2 <: SearchOver end

include("functions.jl")


#cd("C:\\Users\\laidler1\\Documents\\Shapelets\\Code")