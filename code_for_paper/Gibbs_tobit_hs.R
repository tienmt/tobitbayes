library(Rcpp);sourceCpp('test.cpp')
#-------------------------------
# Tobit regression Gibbs sampler
#-------------------------------
gibbs_tobit_horseshoe <- function(y, X, c_sensored = 0, n_iter = 2000, burn_in = 500,a0 = 1, b0 = 1) {
  n <- length(y)
  p <- ncol(X)
  # Initialization
  beta <- rep(0, p)
  sigma2 <- 1
  tau2 <- 1
  xi <- 1
  lambda2 <- rep(1, p)
  nu <- rep(1, p)
  z <- y
  d_i <- y > c_sensored
  n_di <- !d_i 
  total_ndi = sum(n_di)
  tXX <- crossprod(X)
  tX <- t(X)
  # Storage
  beta_store <- matrix(NA, nrow = n_iter, ncol = p)
  hat_sig <- 0
  
  for (iter in 1:n_iter) {
    # Step 1: Sample latent z
    mu <- eigenMapMatMult( X,beta )
    z[n_di] <- rtruncnorm(total_ndi, a = -Inf, b = c_sensored, mean = mu[n_di], sd = sqrt(sigma2))

    # Step 2: Sample beta
    D_inv <- diag(1 / (tau2 * lambda2))
    Sigma_beta <-  chol2inv(chol( tXX / sigma2 + D_inv) )
    tamtam <- tX %*% z
    mu_beta <- eigenMapMatMult( Sigma_beta , tamtam ) / sigma2
    LL <- chol(Sigma_beta)
    beta <- mu_beta + LL %*% rnorm(p)
    
    # Step 3: Sample lambda_j^2
    rate_lambda <- 1 / nu + beta^2 / (2*tau2)
    lambda2 <- sapply(rate_lambda, function(x) 1 / rgamma(1, shape = 1, rate = x ) )
    # Step 4: Sample nu_j
   #  for (j in 1:p) { nu[j] <- 1 / rgamma(1, shape = 1, rate = 1 + 1 / lambda2[j])   }
    nu <- sapply(lambda2, function(x) 1 / rgamma(1, shape = 1, rate = 1 + 1/x )  )
    # Step 5: Sample tau^2
    rate_tau <- 1 / xi + sum(beta^2 / lambda2) / 2
    tau2 <- 1 / rgamma(1, shape = (p + 1)/2, rate = rate_tau)
    
    # Step 6: Sample xi
    xi <- 1 / rgamma(1, shape = 1, rate = 1 + 1 / tau2)
    
    # Step 7: Sample sigma^2
    resid <- z - eigenMapMatMult( X,beta )
    rate_sigma <- b0 + sum(resid^2) / 2
    sigma2 <- 1 / rgamma(1, shape = a0 + n/2, rate = rate_sigma)
    
    if (iter>burn_in) hat_sig <- hat_sig + sigma2/(n_iter-burn_in)
    # Store results
    beta_store[iter, ] <- beta
  }
  # Return posterior samples (after burn-in)
  list(    beta_samples = beta_store[(burn_in + 1):n_iter, ],
            full_beta = beta_store,
            n_eff = n_iter - burn_in,
           sigma2HS = hat_sig)
}
gibbs_tobit_horseshoe <- compiler::cmpfun(gibbs_tobit_horseshoe)




selection_metrics <- function(true_support, selected_support, p = NULL) {
  # true_support: indices of truly relevant variables
  # selected_support: indices selected by your method
  # p: total number of variables (optional, required for TN)
  
  true_support <- as.integer(true_support)
  selected_support <- as.integer(selected_support)
  
  TP <- length(intersect(true_support, selected_support))
  FP <- length(setdiff(selected_support, true_support))
  FN <- length(setdiff(true_support, selected_support))
  if (!is.null(p)) {
    all_indices <- seq_len(p)
    TN <- length(setdiff(all_indices, union(true_support, selected_support)))
  } else {
    TN <- NA  # Cannot compute without knowing p
  }
  
  precision <- if ((TP + FP) == 0) NA else TP / (TP + FP)
  recall    <- if ((TP + FN) == 0) NA else TP / (TP + FN)
  specificity <- if (!is.na(TN) && (TN + FP) > 0) TN / (TN + FP) else NA
  f1 <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) {
    2 * precision * recall / (precision + recall)
  } else {
    NA
  }
  fdr <- if ((TP + FP) == 0) 0 else FP / (TP + FP)
  
  return(c(
    #'Precision' = round( precision, 2 ),
    'TPR' = round( recall,2),
    #'F1' = f1,
    #'Specificity' = specificity,
    'FDR' = round( fdr, 2 )
  ))
}
credible_intervals <- function(samples, level = 0.95) {
  alpha <- (1 - level) / 2
  ci_mat <- apply(samples, 2, function(x) {
    c(Lower = quantile(x, probs = alpha),
      Median = median(x),
      Upper = quantile(x, probs = 1 - alpha))
  })
  ci_mat <- t(ci_mat)
  colnames(ci_mat) <- c("Lower", "Median", "Upper")
  return(ci_mat)
}
