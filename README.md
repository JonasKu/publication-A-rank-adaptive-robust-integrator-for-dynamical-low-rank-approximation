# publication-A-rank-adaptive-robust-integrator-for-dynamical-low-rank-approximation
This code framework can be used to reproduce all numerical results of the paper "A rank-adaptive robust integrator for dynamical low-rank approximation". The code for the radiation transport and uncertainty quantification sections is written in the programming language Julia (Version 1.6.2).

The code uses the following Julia packages: PyPlot, NPZ, ProgressMeter, LinearAlgebra, GSL, FastTransforms, FastGaussQuadrature, DelimitedFiles, LegendrePolynomials, QuadGK, SparseArrays, SphericalHarmonicExpansions,SphericalHarmonics,TypedPolynomials,GSL, MultivariatePolynomials, Einsum. Installing these packages can be done via the Julia package manager (https://docs.julialang.org/en/v1/stdlib/Pkg/).

To run the code, navigate to the corresponding subfolder and open the Julia environment. By typing include("main.jl") the corresponding simulations of the paper are computed and the results will be plotted. 
