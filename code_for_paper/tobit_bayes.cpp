#include <RcppArmadillo.h>
#include <Rmath.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// helper: sample from truncated normal (a, b)
double rtruncnorm(double mean, double sd, double a, double b) {
  double alpha = (a - mean) / sd;
  double beta = (b - mean) / sd;
  double u = R::runif(R::pnorm(alpha, 0.0, 1.0, 1, 0),
                      R::pnorm(beta, 0.0, 1.0, 1, 0));
  double z = R::qnorm(u, 0.0, 1.0, 1, 0);
  return mean + sd * z;
}

// [[Rcpp::export]]
List tobit_bayes_cpp(const arma::vec& y,
                     const arma::mat& X,
                     double c_sensored = 0,
                     int n_iter = 2000,
                     int burn_in = 500,
                     double a0 = 1.0,
                     double b0 = 0.1) {
  
  int n = y.n_elem;
  int p = X.n_cols;
  
  // Initialization
  arma::vec beta = arma::zeros(p);
  double sigma2 = 1.0;
  double tau2 = 1.0;
  double xi = 1.0;
  arma::vec lambda2 = arma::ones(p);
  arma::vec nu = arma::ones(p);
  arma::vec z = y;
  
  uvec d_i = find(y > c_sensored);
  uvec n_di = find(y <= c_sensored);
  int total_ndi = n_di.n_elem;
  
  arma::mat tXX = X.t() * X;
  arma::mat tX = X.t();
  
  arma::mat beta_store = arma::mat(n_iter, p, fill::zeros);
  double hat_sig = 0.1;
  
  arma::vec mu = X * beta;
  
  for (int iter = 0; iter < n_iter; ++iter) {
    
    // --- Step 1: sample latent z
    mu = X * beta;
    for (int i = 0; i < total_ndi; ++i) {
      int idx = n_di[i];
      z[idx] = rtruncnorm(mu[idx], std::sqrt(sigma2), -INFINITY, c_sensored);
    }
    
    // --- Step 2: sample beta using Cholesky solve (faster, more stable)
    arma::mat Q = tXX / sigma2 + arma::diagmat(1.0 / (tau2 * lambda2));  // precision matrix
    arma::mat L = arma::chol(Q, "lower");                                // Q = L * L^T
    
    // Compute mu_beta = Q^{-1} * (tX * z) / sigma2  without forming Q^{-1}
    arma::vec b = (tX * z) / sigma2;
    arma::vec mu_beta = arma::solve(arma::trimatu(L.t()),
                                    arma::solve(arma::trimatl(L), b, arma::solve_opts::fast),
                                    arma::solve_opts::fast);
    
    // Sample from N(mu_beta, Q^{-1})
    arma::vec eps = arma::randn(p);
    arma::vec v = arma::solve(arma::trimatu(L.t()), eps, arma::solve_opts::fast);
    beta = mu_beta + v;
    
    
    // --- Step 3: sample lambda_j^2
    arma::vec beta222 = arma::square(beta);
    arma::vec rate_lambda = 1.0 / nu + beta222 / (2.0 * tau2);
    for (int j = 0; j < p; ++j) {
      lambda2[j] = 1.0 / R::rgamma(1.0, 1.0 / rate_lambda[j]); // shape=1, scale=1/rate

      nu[j] = 1.0 / R::rgamma(1.0, 1.0 / (1.0 + 1.0 / lambda2[j]));
    }
    
    // --- Step 5: sample tau^2
    double rate_tau = 1.0 / xi + arma::sum(beta222 / lambda2) / 2.0;
    tau2 = 1.0 / R::rgamma((p + 1.0) / 2.0, 1.0 / rate_tau);
    
    // --- Step 6: sample xi
    xi = 1.0 / R::rgamma(1.0, 1.0 / (1.0 + 1.0 / tau2));
    
    // --- Step 7: sample sigma^2
    arma::vec resid = z - X * beta;
    double rate_sigma = b0 + arma::dot(resid, resid) / 2.0;
    sigma2 = 1.0 / R::rgamma(a0 + n / 2.0, 1.0 / rate_sigma);
    
    // accumulate posterior mean of sigma²
    if (iter >= burn_in) {
      hat_sig += sigma2 / (n_iter - burn_in);
    }
    
    // store
    beta_store.row(iter) = beta.t();
  }
  
  arma::mat beta_samples = beta_store.rows(burn_in, n_iter - 1);
  
  return List::create(
    Named("beta_samples") = beta_samples,
    Named("full_beta") = beta_store,
    Named("sigma2HS") = hat_sig
  );
}
