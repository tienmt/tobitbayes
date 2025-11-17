# code for reproduce results in the paper
**"High-dimensional Bayesian Tobit regression for censored response with Horseshoe prior."**

Main file for simulation are:
"simu_tobit.R"     for most simulations
"simu_tobit_CIs.R" for simulation assesment fo credible intervals

Other funtions that NEED to be loaded for runing simulations and realdata are:
"tobit_SnS.cpp"         Rcpp function implement of spike and slab prior
"tobit_bayes.cpp"       Rcpp function implement of Horseshoe prior
"tobit_scad_fast.cpp"   Rcpp function implement of SCAD.
"Gibbs_tobit_hs.R"      Auxillary functions

Real data can be re-produced using the code:
"real data tobit.R"
