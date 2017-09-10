# Diffusivity Calculations

These are the fortran and matlab scripts I use to calculate molecular diffusivities (using Tyn and Calus Method).

- To run cd into folder with required input data (e.g. "hephex_experiments") and type in the command line "../source/molecularDiff".  NOTE: this only works with the input file for "hephex_experiments" at the moment, and does not work for the input file for "prpgly_experiments".  The nput file for "progly_experiments" are formated differently.

- The main code is "molecularDiff.Rev02.f90"

- The "liquidDensities.m" script, predicts droplet diameter percent increase due to heating.

- The "decayTimes.m" script calculates viscous decay time.

- The script "molecularDiffusivty_uncertainty_MontiCarlo_xxxx.R" in the folders "hephex_experiments" and "propgly_experiments" calculates D_{AB} and D_{BA} using a Monte Carlo simulation.