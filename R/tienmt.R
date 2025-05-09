#' Gibbs Sampler for Tobit Regression with Horseshoe Prior
#'
#' @param y Response vector (numeric) of length n, possibly censored.
#' @param X Design matrix (numeric), n x p.
#' @param c_sensored Censoring threshold (default 0).
#' @param n_iter Number of MCMC iterations (default 2000).
#' @param burn_in Number of burn-in iterations (default 500).
#' @param a0 Shape parameter for inverse-gamma prior on sigma² (default 1).
#' @param b0 Rate parameter for inverse-gamma prior on sigma² (default 0.1).
#' @param level level for return credible intervals and variable selction, default is 0.95. should be in (0,1).
#'
#' @importFrom stats rnorm rgamma median quantile density acf
#' @importFrom graphics hist lines par
#'
#' @return A list containing:
#' \describe{
#'   \item{beta_samples}{A numeric matrix of dimension \code{(n_iter - burn_in) × p} containing posterior samples of the regression coefficients after burn-in.}
#'   \item{full_beta}{A numeric matrix of dimension \code{n_iter × p} storing all sampled beta values including burn-in iterations. Useful for convergence diagnostics.}
#'   \item{sigma2HS}{A numeric value giving the posterior mean of the error variance \eqn{\sigma^2} over the samples after burn-in.}
#'   \item{selected}{A binary vector giving the selected variable based on credible intervals.}
#' }
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   X <- matrix(rnorm(100*5), 100, 5)
#'   beta0 <- c(2, -1, 0, 0, 0)
#'   y <- X %*% beta0 + rnorm(100)
#'   y[y < 0] <- 0
#'   res <- tobit_bayes(y, X)
#'   (posterior_means <- colMeans(res$beta_samples))
#'   plot(posterior_means, type = 'h', main = "Posterior Means of Beta", ylab = "Mean")
#'   res$selected
#' }
#' @export
tobit_bayes <- function(y, X, c_sensored = 0,
                        n_iter = 2000,
                        burn_in = 500,
                        a0 = 1,
                        b0 = .1,
                        level = 0.95) {
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
  d_i = y > c_sensored
  n_di = !d_i;
  total_ndi = sum(n_di)
  tXX <- crossprod(X)
  tX <- t(X)
  # Storage
  beta_store <- matrix(NA, nrow = n_iter, ncol = p)
  hat_sig <- 0.1

  for (iter in 1:n_iter) {
    # Step 1: Sample latent z
    mu <- X %*%beta
    z[n_di] <- truncnorm::rtruncnorm(total_ndi, a = -Inf, b = c_sensored, mean = mu[n_di], sd = sqrt(sigma2))

    # Step 2: Sample beta
    D_inv <- diag(1 / (tau2 * lambda2))
    Sigma_beta <- chol2inv(chol( tXX / sigma2 + D_inv) )
    tamtam <-  tX %*% z
    mu_beta <-  Sigma_beta %*% tamtam  / sigma2
    #beta <- mvrnorm(1, mu = mu_beta, Sigma = Sigma_beta)
    LL <- chol(Sigma_beta)
    beta <- mu_beta + LL %*% rnorm(p)

    # Step 3: Sample lambda_j^2
    beta222 <- beta^2
    rate_lambda <- 1 / nu + beta222 / (2*tau2)
    lambda2 <- sapply(rate_lambda, function(x) 1 / rgamma(1, shape = 1, rate = x ) )
    # Step 4: Sample nu_j
   #  for (j in 1:p) { nu[j] <- 1 / rgamma(1, shape = 1, rate = 1 + 1 / lambda2[j])   }
    nu <- sapply(lambda2, function(x) 1 / rgamma(1, shape = 1, rate = 1 + 1/x )  )
    # Step 5: Sample tau^2
    rate_tau <- 1 / xi + sum(beta222 / lambda2) / 2
    tau2 <- 1 / rgamma(1, shape = (p + 1)/2, rate = rate_tau)
    # Step 6: Sample xi
    xi <- 1 / stats::rgamma(1, shape = 1, rate = 1 + 1 / tau2)

    # Step 7: Sample sigma^2
    resid <- z - X %*% beta
    rate_sigma <- b0 + sum(resid^2) / 2
    sigma2 <- 1 / rgamma(1, shape = a0 + n/2, rate = rate_sigma)

    if (iter>burn_in) hat_sig <- hat_sig + sigma2/(n_iter-burn_in)
    # Store results
    beta_store[iter, ] <- beta
  }
  credible_intervals <- function(beta_samples, level = 0.95) {
    alpha <- (1 - level) / 2
    ci_mat <- apply(beta_samples, 2, function(x) {
      c(Lower = quantile(x, probs = alpha),
        Median = median(x),
        Upper = quantile(x, probs = 1 - alpha))
    })
    ci_mat <- t(ci_mat)
    colnames(ci_mat) <- c("Lower", "Median", "Upper")
    return(ci_mat)
  }
  ci <- credible_intervals( beta_store[(burn_in + 1):n_iter,]  , level = level)
  selected <- !(ci[, "Lower"] < 0 & ci[, "Upper"] > 0)
  # Return posterior samples (after burn-in)
  list(    beta_samples = beta_store[(burn_in + 1):n_iter, ],
            full_beta = beta_store,
           sigma2HS = hat_sig,
           selected = selected)
}

#' Plot Posterior Distributions of Regression Coefficients
#'
#' @param obj output from (from \code{tobit_bayes}).
#' @param which Optional vector of indices for coefficients to plot.
#' @examples
#' \dontrun{
#'   set.seed(1)
#'   X <- matrix(rnorm(100*5), 100, 5)
#'   beta0 <- c(2, -1, 0, 0, 0)
#'   y <- X %*% beta0 + rnorm(100)
#'   y[y < 0] <- 0
#'   res <- tobit_bayes(y, X)
#'  ( posterior_means <- colMeans(res$beta_samples))
#'   plot_posterior(res)
#' }
#' @export
plot_posterior <- function(obj , which = 1:ncol(beta_samples)) {
  beta_samples <- obj$beta_samples
  op <- par(mfrow = c(ceiling(length(which) / 2), 2))
  on.exit(par(op))
  for (j in which) {
    hist(beta_samples[, j], main = paste("Beta", j), xlab = "", col = "lightblue", border = "white", probability = TRUE)
    lines(density(beta_samples[, j]), col = "darkblue", lwd = 2)
  }
}

#' Posterior Predictive Means from Tobit Model
#'
#' Computes the posterior predictive mean for new observations using samples of regression coefficients.
#'
#' @param beta_samples A matrix of posterior samples for regression coefficients (e.g., from \code{tobit_bayes}).
#' @param Xnew A matrix of new predictor values. Must have the same number of columns as the original design matrix.
#'
#' @return A numeric vector of posterior predictive means for each row of \code{Xnew}.
#'
#' @examples
#' \dontrun{
#'   set.seed(42)
#'   X <- matrix(rnorm(100 * 3), 100, 3)
#'   beta0 <- c(1, 0, 2)
#'   y <- X %*% beta0 + rnorm(100)
#'   y[y < 0] <- 0
#'   fit <- tobit_bayes(y, X)
#'   Xnew <- matrix(rnorm(5 * 3), 5, 3)
#'   preds <- predict_tobit_hs(fit$beta_samples, Xnew)
#'   print(preds)
#' }
#' @export
predict_tobit_hs <- function(beta_samples, Xnew) {
  apply(beta_samples, 1, function(beta) Xnew %*% beta) |> rowMeans()
}

#' Effective Sample Size (ESS) Estimation
#'
#' Approximates the effective sample size for each regression coefficient based on autocorrelation.
#'
#' @param beta_samples A matrix of posterior samples (from \code{tobit_bayes}).
#'
#' @return A numeric vector giving the effective sample size for each coefficient.
#'
#' @examples
#' \dontrun{
#'   set.seed(123)
#'   X <- matrix(rnorm(200), 100, 2)
#'   beta0 <- c(1, 0)
#'   y <- X %*% beta0 + rnorm(100)
#'   y[y < 0] <- 0
#'   res <- tobit_bayes(y, X)
#'   ess <- effective_sample_size(res$beta_samples)
#'   print(ess)
#' }
#' @export
effective_sample_size <- function(beta_samples) {
  apply(beta_samples, 2, function(x) {
    n <- length(x)
    acf_vals <- acf(x, lag.max = min(100, n - 1), plot = FALSE)$acf[-1]
    rho <- acf_vals[which(acf_vals > 0)]
    ess <- n / (1 + 2 * sum(rho))
    round(ess)
  })
}




