#include <RcppArmadillo.h>
#include <Rmath.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// helper: truncated normal sampler
double rtruncnorm(double mean, double sd, double a, double b) {
  double alpha = (a - mean) / sd;
  double beta = (b - mean) / sd;
  double u = R::runif(R::pnorm(alpha, 0.0, 1.0, 1, 0),
                      R::pnorm(beta, 0.0, 1.0, 1, 0));
  double z = R::qnorm(u, 0.0, 1.0, 1, 0);
  return mean + sd * z;
}

// [[Rcpp::export]]
List tobit_spikeslab_cpp(const arma::vec& y,
                         const arma::mat& X,
                         double c_sensored = 0.0,
                         int n_iter = 2000,
                         int burn_in = 500,
                         double a0 = 1.0,
                         double b0 = 0.1,
                         double tau0 = 0.01,     // spike variance
                         double tau1 = 1.0,      // slab variance
                         double a_pi = 1.0,      // Beta prior for pi
                         double b_pi = 1.0) {
  
  int n = y.n_elem;
  int p = X.n_cols;
  
  // --- Initialization ---
  arma::vec beta = arma::zeros(p);
  double sigma2 = 1.0;
  arma::vec z = y;
  arma::vec gamma = arma::ones(p);  // inclusion indicators
  double pi = 0.5;                  // prior inclusion probability
  
  uvec d_i = find(y > c_sensored);
  uvec n_di = find(y <= c_sensored);
  int total_ndi = n_di.n_elem;
  
  arma::mat tXX = X.t() * X;
  arma::mat tX = X.t();
  
  arma::mat beta_store(n_iter, p, fill::zeros);
  arma::mat gamma_store(n_iter, p, fill::zeros);
  arma::vec sigma_store(n_iter, fill::zeros);
  arma::vec pi_store(n_iter, fill::zeros);
  
  arma::vec mu = X * beta;
  
  // Gibbs sampling
  for (int iter = 0; iter < n_iter; ++iter) {
    
    // --- Step 1: sample latent z_i
    mu = X * beta;
    for (int i = 0; i < total_ndi; ++i) {
      int idx = n_di[i];
      z[idx] = rtruncnorm(mu[idx], std::sqrt(sigma2), -INFINITY, c_sensored);
    }
    
    // --- Step 2: sample beta | rest using Cholesky solve (no inverse)
    
    // Construct diagonal prior precision D_inv efficiently
    arma::vec Dinv_diag(p);
    for (int j = 0; j < p; ++j) {
      double tau_j2 = (gamma[j] == 1.0) ? (tau1 * tau1) : (tau0 * tau0);
      Dinv_diag[j] = 1.0 / (tau_j2 * sigma2);
    }
    
    // Precision matrix Q = X'X / sigma2 + D_inv
    arma::mat Q = tXX / sigma2 + arma::diagmat(Dinv_diag);
    
    // Cholesky factorization: Q = L * L^T
    arma::mat L = arma::chol(Q, "lower");
    
    // Compute mu_beta = Q^{-1} * (X' z / sigma2) without forming inverse
    arma::vec b = (tX * z) / sigma2;
    arma::vec mu_beta = arma::solve(arma::trimatu(L.t()),
                                    arma::solve(arma::trimatl(L), b, arma::solve_opts::fast),
                                    arma::solve_opts::fast);
    
    // Sample beta ~ N(mu_beta, Q^{-1}) via L^{-T} * eps
    arma::vec eps = arma::randn(p);
    arma::vec v = arma::solve(arma::trimatu(L.t()), eps, arma::solve_opts::fast);
    beta = mu_beta + v;
    
    
    
    
    
    // --- Step 3: sample gamma_j | beta_j, pi
    for (int j = 0; j < p; ++j) {
      double dens1 = R::dnorm(beta[j], 0.0, std::sqrt(sigma2 * tau1 * tau1), 0);
      double dens0 = R::dnorm(beta[j], 0.0, std::sqrt(sigma2 * tau0 * tau0), 0);
      double p1 = pi * dens1;
      double p0 = (1.0 - pi) * dens0;
      double prob = p1 / (p1 + p0 + 1e-12);
      gamma[j] = R::rbinom(1, prob);
    }
    
    // --- Step 4: sample pi | gamma
    double a_post = a_pi + arma::sum(gamma);
    double b_post = b_pi + p - arma::sum(gamma);
    pi = R::rbeta(a_post, b_post);
    
    // --- Step 5: sample sigma2 | z, beta
    arma::vec resid = z - X * beta;
    double shape_sigma = a0 + n / 2.0;
    double rate_sigma = b0 + arma::dot(resid, resid) / 2.0;
    sigma2 = 1.0 / R::rgamma(shape_sigma, 1.0 / rate_sigma);
    
    // --- Store samples
    beta_store.row(iter) = beta.t();
    gamma_store.row(iter) = gamma.t();
    sigma_store[iter] = sigma2;
    pi_store[iter] = pi;
  }
  
  arma::mat beta_samples = beta_store.rows(burn_in, n_iter - 1);
  arma::mat gamma_samples = gamma_store.rows(burn_in, n_iter - 1);
  
  return List::create(
    Named("beta_samples") = beta_samples,
    Named("gamma_samples") = gamma_samples,
    Named("sigma2_samples") = sigma_store.rows(burn_in, n_iter - 1),
    Named("pi_samples") = pi_store.rows(burn_in, n_iter - 1),
    Named("mean_beta") = mean(beta_samples, 0),
    Named("inclusion_prob") = mean(gamma_samples, 0)
  );
}
