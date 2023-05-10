# Simulation Shapelets

[![DOI](https://zenodo.org/badge/631010467.svg)](https://zenodo.org/badge/latestdoi/631010467)

This repository contains data and code used for a shapelet analysis of simulation trajectories.

To use, clone the repository. Refer to the [instructions](https://pkgdocs.julialang.org/v1.2/environments/#Using-someone-else's-project-1) for setting up Julia environments.

The following files should be used to run the analysis:
- Breakdowns.jl: Perform analysis to identify the impact of machine breakdowns on the dynamic throughput of a wafer fab simulation.
- Validation.jl: Perform analysis for dynamic model validation of two models of a tandem queue.
- Multivariate.jl: Perform analysis on the multivariate system state of a tandem queue.
- figures.jl: Generate figures from the results of the previous three files.
