% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tienmt.R
\name{plot_posterior}
\alias{plot_posterior}
\title{Plot Posterior Distributions of Regression Coefficients}
\usage{
plot_posterior(obj, which = 1:ncol(beta_samples))
}
\arguments{
\item{obj}{output from (from \code{tobit_bayes}).}

\item{which}{Optional vector of indices for coefficients to plot.}
}
\description{
Plot Posterior Distributions of Regression Coefficients
}
\examples{
\dontrun{
  set.seed(1)
  X <- matrix(rnorm(100*5), 100, 5)
  beta0 <- c(2, -1, 0, 0, 0)
  y <- X \%*\% beta0 + rnorm(100)
  y[y < 0] <- 0
  res <- tobit_bayes(y, X)
 ( posterior_means <- colMeans(res$beta_samples))
  plot_posterior(res)
}
}
