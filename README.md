# tobitbayes
Bayesian sparse Tobit regression for censored response via Gibbs sampler. This is based the paper: "High-dimensional Bayesian Tobit regression for censored response with Horseshoe prior".

Install by using:

devtools::install_github('tienmt/tobitbayes')

Example:

library(tobitbayes)
 set.seed(1)
   X <- matrix(rnorm(100*5), 100, 5)
   beta0 <- c(2, -1, 0, 0, 0)
   y <- X %*% beta0 + rnorm(100)
   y[y < 0] <- 0
   res <- tobit_bayes(y, X)
   (posterior_means <- colMeans(res$beta_samples))
   plot(posterior_means, type = 'h', main = "Posterior Means of Beta", ylab = "Mean")
  res$selected
