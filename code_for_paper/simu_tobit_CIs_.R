###############################################################
#         Simulation + Gibbs Sampler for Tobit Regression
#      Comparing Horseshoe (HS) vs Spike-and-Slab (SnS) Priors
#                 Fully commented and reader-friendly
###############################################################

# -------------------------------------------------------------
# Load required libraries
# -------------------------------------------------------------
# devtools::install_github("TateJacobson/tobitnet")  # Optional dependency

library(tobitnet)   # Tobit regression utilities (if needed)
library(tictoc)     # Simple timing utility
library(truncnorm)  # Sampling from truncated normal distributions
library(MASS)       # Multivariate normal generator (mvrnorm)
library(Rcpp)       # Needed for calling C++-compiled samplers

# Load Rcpp implementations of the samplers
sourceCpp('tobit_bayes.cpp')         # Horseshoe prior Gibbs sampler (C++)
sourceCpp('tobit_scad_fast.cpp')     # SCAD implementation (optional)
sourceCpp('tobit_SnS.cpp')           # Spike-and-slab Gibbs sampler (C++)

# For comparison with your R Gibbs sampler (not used here)
source('Gibbs_tobit_hs.R')

# -------------------------------------------------------------
# Helper: credible intervals for posterior samples
# -------------------------------------------------------------
# Computes median + (Lower, Upper) credible interval for each parameter.
credible_intervals <- function(samples, level = 0.95) {
  alpha <- (1 - level) / 2
  
  # Apply over columns (each column = MCMC samples for one coefficient)
  ci_mat <- apply(samples, 2, function(x) {
    c(
      Lower  = quantile(x, probs = alpha),
      Median = median(x),
      Upper  = quantile(x, probs = 1 - alpha)
    )
  })
  
  ci_mat <- t(ci_mat)  # Make each row a coefficient
  colnames(ci_mat) <- c("Lower", "Median", "Upper")
  return(ci_mat)
}


###############################################################
#                      Simulation Settings
###############################################################

n <- 100            # Training sample size
n_test <- 1000      # Test sample size
n_all <- n + n_test # Total data generated

p <- 120            # Number of predictors
s0 <- 4             # Number of non-zero coefficients

# Truth: first s0/2 = 1, last s0/2 = -1
beta_true <- rep(0, p)
beta_true[1:s0] <- c(rep(1, s0/2), rep(-1, s0/2))

# Covariance structure with correlation rho^{|i-j|}
rho <- 0
Sigma <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
LL <- chol(Sigma)   # Cholesky for generating correlated predictors


###############################################################
#              Storage for Monte Carlo replications
###############################################################

HSout <- vector("list", 100)     # Horseshoe results
SnS_out <- vector("list", 100)   # Spike-slab results


###############################################################
#             Monte Carlo Simulation Loop (100 reps)
###############################################################

for (ss in 1:100) {
  
  # -----------------------------------------------------------
  # Generate design matrix X (training) and X_test (testing)
  # -----------------------------------------------------------
  xall <- matrix(rnorm(n_all * p), n_all, p) %*% LL
  X     <- xall[1:n, ]
  Xtest <- xall[-(1:n), ]
  
  # -----------------------------------------------------------
  # Generate latent Tobit response:  z = Xβ + ε
  # Tobit observes y = max(censored, z)
  # -----------------------------------------------------------
  z_true <- xall %*% beta_true + rnorm(n_all, sd = 1)
  
  censored <- 0   # Left-censoring threshold
  y     <- pmax(censored, z_true[1:n])
  ytest <- pmax(censored, z_true[-(1:n)])
  
  # -----------------------------------------------------------
  # 1. Run Horseshoe Gibbs sampler (C++)
  # -----------------------------------------------------------
  result_cpp <- tobit_bayes_cpp(
    y, X,
    c_sensored = censored,
    n_iter     = 5000,
    burn_in    = 1000
  )
  
  # Posterior CIs
  ci_hs <- credible_intervals(result_cpp$beta_samples)
  
  CIhs_lower  <- ci_hs[, "Lower"]
  CIhs_upper  <- ci_hs[, "Upper"]
  CIhs_length <- CIhs_upper - CIhs_lower
  
  # Store HS: (average CI length, coverage rate)
  HSout[[ss]] <- c(
    mean(CIhs_length),
    mean(beta_true >= CIhs_lower & beta_true <= CIhs_upper)
  )
  
  # -----------------------------------------------------------
  # 2. Run Spike-and-Slab Gibbs sampler (C++)
  # -----------------------------------------------------------
  res_SnS_cpp <- tobit_spikeslab_cpp(
    y, X,
    c_sensored = censored,
    n_iter     = 5000,
    burn_in    = 1000
  )
  
  ci_Sns <- credible_intervals(res_SnS_cpp$beta_samples)
  
  CIss_lower  <- ci_Sns[, "Lower"]
  CIss_upper  <- ci_Sns[, "Upper"]
  CIss_length <- CIss_upper - CIss_lower
  
  # FIXED: previously you stored mean(CIhs_length) by accident
  SnS_out[[ss]] <- c(
    mean(CIss_length),
    mean(beta_true >= CIss_lower & beta_true <= CIss_upper)
  )
  # Progress
  print(ss)
}

# Clean-up
rm(Sigma, LL, X, xall)

# Save full workspace for later analysis
save.image(file = "tobit_p300n200_s32_.rda")
