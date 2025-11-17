

# Code for Reproducing the Results in

**“High-Dimensional Bayesian Tobit Regression for Censored Response with the Horseshoe Prior.”**

This repository contains all scripts and functions required to reproduce the simulation studies and real-data analysis presented in the paper.

## Main Simulation Scripts

* **`simu_tobit.R`** – Primary script for running most simulation experiments.
* **`simu_tobit_CIs.R`** – Script for evaluating the performance of credible intervals.

## Required Supporting Functions

The following files must be loaded before running the simulations or real-data analysis:

* **`tobit_SnS.cpp`** – Rcpp implementation of the spike-and-slab prior.
* **`tobit_bayes.cpp`** – Rcpp implementation of the Horseshoe prior.
* **`tobit_scad_fast.cpp`** – Rcpp implementation of the SCAD penalty.
* **`Gibbs_tobit_hs.R`** – Additional helper functions for the Gibbs sampler.

## Real Data Analysis

The real-data results in the paper can be reproduced using:

* **`real data tobit.R`**


