# tobitbayes

Bayesian sparse Tobit regression for censored response via Gibbs sampler. This is based on the paper:  
**"High-dimensional Bayesian Tobit regression for censored response with Horseshoe prior."**

## Installation

Install the package using:

```r
devtools::install_github('tienmt/tobitbayes')

library(tobitbayes)

# simulate data
set.seed(1)
X <- matrix(rnorm(100 * 5), 100, 5)
beta0 <- c(2, -1, 0, 0, 0)
y <- X %*% beta0 + rnorm(100)

# censor the response to value c = 0
y[y < 0] <- 0

# fit the Gibbs sampler
res <- tobit_bayes(y, X)

# get the posterior mean and compare it to the true beta0
(posterior_means <- colMeans(res$beta_samples))
# plot output
plot(posterior_means, type = 'h',
     main = "Posterior Means of Beta",
     ylab = "Mean")

### get the selected variables from the Gibbs sampler
res$selected
