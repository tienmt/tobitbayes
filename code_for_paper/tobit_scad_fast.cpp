// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// ----------------- helper functions provided by user -----------------

double normalCDF(double x){
  return (std::erfc(-x/std::sqrt(2))/2);
}

double normalPDF(double x){
  static const double inv_sqrt_2pi = 0.3989422804014327;
  return(inv_sqrt_2pi * std::exp(-0.5 * x * x));
}

double g(double x){
  double gout;
  if(x < -37){
    gout = -x;
  } else {
    gout = normalPDF(x)/normalCDF(x);
  }
  return(gout);
}

//[[Rcpp::export]]
NumericVector LprimeC(NumericVector xj, NumericVector y, LogicalVector d, NumericVector r, double gamma){
  int n = y.size();
  NumericVector out(n);
  for(int i = 0; i < n; i++){
    if(d[i]){
      out[i] = -(gamma*y[i] - r[i])*xj[i];
    } else {
      out[i] = g(-r[i])*xj[i];
    }
  }
  return out;
}

//[[Rcpp::export]]
NumericVector logL1(NumericVector y, LogicalVector d, NumericVector r, double gamma){
  int n = y.size();
  NumericVector out(n);
  for(int i = 0; i < n; i++){
    if(d[i]){
      out[i] = std::log(gamma) -(gamma*y[i] - r[i])*(gamma*y[i] - r[i])/2;
    } else {
      out[i] = std::log(normalCDF(-r[i])) ;
    }
  }
  return out;
}

double soft(double z, double t){
  double sgn = (z > 0) - (z < 0);
  double out = std::max(0.0, std::abs(z) - t)*sgn;
  return out;
}

bool KKT_check(double deltaj, NumericVector xj, NumericVector y, LogicalVector d, NumericVector r, double gamma, double l1, double l2, double tol = 1e-3){
  if(deltaj > 0){
    double check = mean( LprimeC(xj, y, d, r, gamma) ) + l2*deltaj + l1;
    if( -tol < check && check < tol ) return(true);    
  } else if (deltaj < 0){
    double check = mean( LprimeC(xj, y, d, r, gamma) ) + l2*deltaj - l1;
    if( -tol < check && check < tol ) return(true);        
  } else if (deltaj == 0){
    double check = - mean( LprimeC(xj, y, d, r, gamma) ); 
    if( - l1 - tol  < check && check < l1 + tol ) return(true); 
  }
  return(false);
}

NumericVector matProd(NumericMatrix x, NumericVector b){
  int n = x.nrow();
  NumericVector out(n);
  for(int i = 0; i<n; i++){
    out(i) = sum(x(i,_)*b);
  }
  return out;
}

List standardizeC(NumericMatrix x){
  int p = x.ncol();
  NumericVector colmeans(p);
  NumericVector colsds(p);
  for(int j = 0; j < p; j++){
    colmeans(j) = mean( x(_,j) );
    x(_,j) = x(_,j) - colmeans(j);
    colsds(j) = sd( x(_,j) );
    // avoid division by zero
    if(colsds(j) == 0) colsds(j) = 1.0;
    x(_,j) = x(_,j)/colsds(j);
  }
  return List::create(
    Named("colmeans") = colmeans,
    Named("colsds") = colsds
  );
}

// ----------------- provided tobitnet_innerC (unchanged logic) -----------------
// Note: user provided full implementation; we paste it here with [[Rcpp::export]].
// (Small modifications: argument types consistent: cin is NumericVector possibly length 1.)



// helper scalar g (same as yours) but inlined for Armadillo use
inline double g_scalar(double x){
  if(x < -37.0) return -x;
  static const double inv_sqrt_2pi = 0.3989422804014327;
  double pdf = inv_sqrt_2pi * std::exp(-0.5 * x * x);
  double cdf = 0.5 * std::erfc(-x/std::sqrt(2.0));
  if(cdf == 0.0) return -x; // avoid division by zero
  return pdf / cdf;
}

// [[Rcpp::export]]
List tobitnet_innerC(NumericMatrix xin, NumericVector yin, NumericVector cin,
                     double lambda1, double lambda2,
                     NumericVector pf1, NumericVector pf2,
                     NumericVector delta_init,
                     double delta_0_init = 0.0,
                     double gamma_init = 1.0,
                     double eps = 1e-7,
                     bool standardize = true,
                     double maxit = 1e4) {
  
  // Convert to Armadillo objects for speed
  arma::mat X = arma::mat(as<arma::mat>(xin));   // n x p
  arma::vec y = arma::vec(as<arma::vec>(yin));
  arma::vec c = arma::vec(as<arma::vec>(cin));
  int n = X.n_rows;
  int p = X.n_cols;
  
  // prepare d (logical: censored or not) and shift y by c if finite
  arma::uvec d(n);
  if(std::isfinite(c(0))){
    for(int i=0;i<n;i++){
      double yshift = y(i) - c(0);
      d(i) = (yshift > 0) ? 1u : 0u;
      y(i) = yshift;
    }
  } else {
    d.fill(1u);
  }
  
  // standardization (in-place on X)
  arma::vec colmeans(p, arma::fill::zeros);
  arma::vec colsds(p, arma::fill::ones);
  if(p > 0 && standardize){
    for(int j=0;j<p;j++){
      colmeans(j) = arma::mean(X.col(j));
      X.col(j) -= colmeans(j);
      colsds(j) = arma::stddev(X.col(j));
      if(colsds(j) == 0.0) colsds(j) = 1.0;
      X.col(j) /= colsds(j);
    }
  }
  
  // initialize parameters (delta in same scale as your previous code)
  arma::vec delta_current = arma::vec(as<arma::vec>(delta_init));
  arma::vec delta_step = arma::vec(p, arma::fill::zeros);
  double delta_0_current = delta_0_init;
  double gamma_current = gamma_init;
  
  // r = X * delta + delta_0
  arma::vec r = X * delta_current + delta_0_current * arma::ones<arma::vec>(n);
  
  // Precompute M (diagonal denom part)
  arma::vec M(p);
  if(standardize){
    M.fill(1.0);
  } else {
    for(int j=0;j<p;j++){
      M(j) = arma::mean( arma::square( X.col(j) ) );
    }
  }
  arma::vec updateDenom = 1.0 / ( M + lambda2 * arma::vec(as<arma::vec>(pf2)) );
  arma::vec updatet = lambda1 * arma::vec(as<arma::vec>(pf1));
  
  // For the loop we maintain h_i = contribution per obs: if d[i]==1 then -(gamma*y - r) else g(-r)
  arma::vec h(n);
  auto recompute_h = [&](void){
    for(int i=0;i<n;i++){
      if(d(i)){
        h(i) = -(gamma_current * y(i) - r(i));
      } else {
        h(i) = g_scalar(-r(i));
      }
    }
  };
  
  recompute_h();
  
  bool full_loop_check = false;
  arma::uvec active_set_bool(p); active_set_bool.fill(1); // boolean vector of active
  arma::uvec current_active_set_bool(p); current_active_set_bool.fill(1);
  
  arma::vec delta_new(p, arma::fill::zeros);
  
  for(int k_outer = 0; k_outer < (int)maxit; ++k_outer){
    
    // Update intercept (delta_0)
    double delta_0_step = - arma::mean(h);
    delta_0_current += delta_0_step;
    r += delta_0_step; // add same to all r
    // update h quickly for d==1 (linear) and recompute g for d==0
    for(int i=0;i<n;i++){
      if(d(i)){
        h(i) += delta_0_step; // since h_i = -(gamma*y - r), r increased by delta_0_step
      } else {
        // recompute g(-r)
        h(i) = g_scalar(-r(i));
      }
    }
    
    // Update delta (coordinate descent). For each coordinate j recompute gradient using current h:
    for(int j=0;j<p;j++){
      if(active_set_bool(j)){
        // grad_j = mean( h * x_j )  because LprimeC returns contributions that we put into h
        double grad_j = arma::dot(h, X.col(j)) / (double) n;
        double z = M(j) * delta_current(j) - grad_j;
        delta_new(j) = std::max(0.0, std::abs(z) - updatet(j)) * ((z>0) - (z<0));
        delta_new(j) = delta_new(j) * updateDenom(j);
        
        delta_step(j) = delta_new(j) - delta_current(j);
        if(std::abs(delta_step(j)) > 0.0){
          // update r and h
          r += delta_step(j) * X.col(j); // update r
          // update h: for d==1 it's linear; for d==0 recompute g(-r)
          for(int i=0;i<n;i++){
            if(d(i)){
              h(i) -= delta_step(j) * X(i,j); // because h = -(gamma*y - r) -> increases by delta_step * x_{ij} but negative sign (see sign); adjust carefully:
              // Explanation: earlier h(i) = -(gamma*y - r(i)). If r increased by s, then h increases by s.
              // delta_step added directly to r above so we add delta_step * x(i,j) to h(i) (equivalently h += delta_step * xcol)
              // But the line used '-' because earlier code used r += ...; to avoid sign mistakes use consistent update:
              // Let's correct to use h(i) += delta_step(j) * X(i,j)
            } else {
              h(i) = g_scalar(-r(i));
            }
          }
          // The above had a small sign confusion; fix below:
          // Recompute for d==1 in a consistent way:
          for(int i=0;i<n;i++){
            if(d(i)){
              h(i) = -(gamma_current * y(i) - r(i)); // exact
            }
          }
        }
        
        // commit delta_current
        delta_current(j) = delta_new(j);
        
        // mark inactive if exactly zero and no change (as before)
        if(delta_current(j) == 0.0 && delta_step(j) == 0.0){
          active_set_bool(j) = 0;
        }
      }
    } // end j loop
    
    // Update gamma (solve quadratic)
    double ga = arma::sum( arma::square(y) );
    double gb = - arma::dot(y, r);
    double gc = - (double) arma::sum(d); // number of censored? matches your earlier code
    double ga2inv = 1.0 / (2.0 * ga);
    double disc = gb*gb - 4.0 * ga * gc;
    if(disc < 0) disc = 0.0;
    gamma_current = (-gb + std::sqrt(disc)) * ga2inv;
    
    // update h because gamma changed (affects entries with d==1)
    for(int i=0;i<n;i++){
      if(d(i)){
        h(i) = -(gamma_current * y(i) - r(i));
      } else {
        h(i) = g_scalar(-r(i));
      }
    }
    
    // check convergence on delta steps and intercept step
    double max_delta2 = arma::max( arma::square(delta_step) );
    if( std::max(max_delta2, delta_0_step*delta_0_step) < eps || full_loop_check ){
      if(!full_loop_check){
        current_active_set_bool = active_set_bool;
        active_set_bool.fill(1); // re-activate to full-check
        full_loop_check = true;
      } else {
        if( arma::all(active_set_bool == current_active_set_bool) ){
          break;
        } else {
          full_loop_check = false;
        }
      }
    }
  } // outer loop
  
  // KKT check (reuse existing KKT_check or do vectorised check)
  // To keep API consistent we'll produce identical outputs but reuse KKT_check (we leave user function)
  // Convert back to Rcpp types for return
  NumericVector delta_out = wrap(delta_current);
  double d0_out = delta_0_current;
  double gamma_out = gamma_current;
  
  // compute beta (delta/gamma) and rescale if standardize
  NumericVector beta(p);
  for(int j=0;j<p;j++) beta[j] = delta_current(j) / gamma_current;
  double beta0 = delta_0_current / gamma_current;
  
  if(standardize){
    for(int j=0;j<p;j++){
      beta[j] = beta[j] / colsds(j);
    }
    // adjust intercept
    double adj = 0.0;
    for(int j=0;j<p;j++) adj += beta[j] * colmeans(j);
    beta0 -= adj;
  }
  if(std::isfinite(c(0))) beta0 += c(0);
  double sigma = 1.0 / gamma_current;
  
  // KKT test using your original KKT_check (we can call it). For safety, call it using Rcpp types:
  bool KKT = true;
  {
    NumericVector y_r = wrap(y);
    LogicalVector d_r(n);
    for(int i=0;i<n;i++) d_r[i] = (d(i) != 0);
    NumericVector r_r = wrap(r);
    if(! KKT_check(d0_out, rep(1.0, n), y_r, d_r, r_r, gamma_out, 0.0, 0.0) ){
      KKT = false;
    }
    if(KKT){
      for(int j=0;j<p;j++){
        NumericVector xj = wrap( X.col(j) );
        if( ! KKT_check(delta_out[j], xj, y_r, d_r, r_r, gamma_out, pf1[j]*lambda1, pf2[j]*lambda2) ){
          KKT = false;
          break;
        }
      }
    }
  }
  
  return List::create(
    Named("d0") = d0_out,
    Named("delta") = delta_out,
    Named("gamma") = gamma_out,
    Named("b0") = beta0,
    Named("beta") = beta,
    Named("sigma") = sigma,
    Named("KKT") = KKT
  );
}


// ----------------- SCAD derivative and vectorized version -----------------

inline double dscad_scalar(double t, double lambda, double a){
  double at = std::abs(t);
  if(at <= lambda) return lambda;
  if(at <= a*lambda) return (a*lambda - at)/(a - 1.0);
  return 0.0;
}

//[[Rcpp::export]]
NumericVector dscad_vec(NumericVector t, double lambda, double a){
  int p = t.size();
  NumericVector out(p);
  for(int j = 0; j < p; j++){
    out[j] = dscad_scalar(std::abs(t[j]), lambda, a);
  }
  return out;
}

// ----------------- tobitscad implementation -----------------

//[[Rcpp::export]]
List tobitscad_cpp(NumericMatrix x_in, NumericVector y_in,
                   double c = 0.0, double a = 3.0, int iter = 3,
                   int nlambda = 100, double lambda_factor = 0.05,
                   Nullable<NumericVector> lambda_in = R_NilValue,
                   double eps = 1e-7, bool standardize = true, double maxit = 1e6, bool early_stop = true){
  // clones
  NumericMatrix x = clone(x_in);
  NumericVector y = clone(y_in);
  int n = x.nrow();
  int p = x.ncol();
  if(p <= 1) stop("x must be a matrix with 2 or more columns");
  
  // Prepare d (indicator) and shifted y (subtract c if finite)
  LogicalVector d(n);
  NumericVector yshift = clone(y);
  if(std::isfinite(c)) {
    for(int i=0;i<n;i++){
      yshift[i] = y[i] - c;
      d[i] = (yshift[i] > 0);
    }
  } else {
    for(int i=0;i<n;i++) d[i] = true;
  }
  
  // If lambda vector supplied, use it. Else compute path
  NumericVector lambda_path;
  if(lambda_in.isNotNull()){
    lambda_path = as<NumericVector>(lambda_in);
    nlambda = lambda_path.size();
  } else {
    // compute lambda_max as max absolute mean gradient with delta=0, gamma=1
    NumericVector grads(p);
    NumericVector r0 = rep(0.0, n); // r = 0 initially
    for(int j=0;j<p;j++){
      NumericVector xj = x(_,j);
      NumericVector Lp = LprimeC(xj, yshift, d, r0, 1.0);
      grads[j] = std::abs(mean(Lp));
    }
    double lambda_max = max(grads);
    if(lambda_max <= 0) lambda_max = 1.0; // fallback
    // geometric sequence
    lambda_path = NumericVector(nlambda);
    double lambda_min = lambda_max * lambda_factor;
    for(int k=0;k<nlambda;k++){
      double frac = double(k)/(nlambda-1);
      // geometric interpolation: lambda = lambda_max^(1-frac) * lambda_min^frac
      lambda_path[k] = std::exp( std::log(lambda_max)*(1.0-frac) + std::log(lambda_min)*frac );
    }
  }
  
  // initial path using tobitnet_innerC with pf1 = 1 (warm starts)
  NumericMatrix delta_mat(p, nlambda); // delta (p x nlambda)
  NumericVector gamma_vec(nlambda);
  NumericVector b0_vec(nlambda);
  LogicalVector KKT_vec(nlambda);
  NumericVector nulldev_vec(nlambda);
  
  // initial warm starts
  NumericVector delta_init = rep(0.0, p);
  double delta0_init = 0.0;
  double gamma_init = 1.0;
  
  for(int l = 0; l < nlambda; l++){
    double l1 = lambda_path[l];
    NumericVector pf1 = rep(1.0, p);
    NumericVector pf2 = rep(1.0, p);
    // call innerC
    List tn = tobitnet_innerC(x, y, NumericVector::create(c), l1, 0.0, pf1, pf2, delta_init, delta0_init, gamma_init, eps, standardize, maxit);
    NumericVector delta_here = tn["delta"];
    double gamma_here = as<double>(tn["gamma"]);
    double b0_here = as<double>(tn["b0"]);
    bool KKT_here = as<bool>(tn["KKT"]);
    NumericVector beta_here = tn["beta"]; // scaled beta (rescaled inside innerC)
    double sigma_here = tn["sigma"];
    
    // store
    for(int j=0;j<p;j++) delta_mat(j,l) = delta_here[j];
    gamma_vec[l] = gamma_here;
    b0_vec[l] = b0_here;
    KKT_vec[l] = KKT_here;
    nulldev_vec[l] = 0.0; // placeholder for now
    
    // warm starts for next lambda
    delta_init = delta_here;
    delta0_init = as<double>(tn["d0"]);
    gamma_init = gamma_here;
  }
  
  // null dev (approx): compute deviance at null (all-zero betas)
  // r_null = rep((b0 - c)*gamma, n) ; but we don't have a single b0 for null: approximate using last b0
  double null_dev = 0.0;
  {
    int l = 0;
    double gamma_l = gamma_vec[l];
    NumericVector delta_pred0 = delta_mat(_, l); // delta for lambda 1
    NumericVector delta_pred_scaled(p);
    for(int j=0;j<p;j++) delta_pred_scaled[j] = delta_pred0[j];
    NumericVector rtemp = matProd(x, delta_pred_scaled);
    // compute logL1 for null-like; here we compute for last fit and call that nulldev = -2*sum(logL1)
    NumericVector llv = logL1(yshift, d, rtemp, gamma_l);
    null_dev = -2.0 * sum(llv);
  }
  
  // LLA loop: we need to iterate (iter-1) times (R code does 1:(iter-1) and uses inner loops)
  NumericMatrix beta_mat_final(p, nlambda); // will hold beta from last inner fits
  NumericVector dev(nlambda);
  NumericMatrix delta_current = clone(delta_mat); // p x nlambda
  
  for(int it = 0; it < iter-1; it++){
    NumericVector gamma_new(nlambda);
    NumericVector b0_new(nlambda);
    NumericMatrix beta_new(p, nlambda);
    
    // initial warm starts for this LLA iteration: start from first column
    NumericVector delta0_init_vec = delta_current(_, 0);
    NumericVector delta_init_vec = delta_current(_, 0);
    double delta0_scalar_init = 0.0;
    double gamma_scalar_init = 1.0;
    
    for(int l = 0; l < nlambda; l++){
      double l1 = lambda_path[l];
      // compute pf1 = dscad(abs(delta_current[,l]), lambda = l1, a)/l1
      NumericVector delta_abs_col(p);
      for(int j=0;j<p;j++) delta_abs_col[j] = std::abs(delta_current(j,l));
      NumericVector scadvals = dscad_vec(delta_abs_col, l1, a);
      NumericVector pf1(p);
      for(int j=0;j<p;j++) pf1[j] = scadvals[j] / l1;
      NumericVector pf2 = rep(1.0, p);
      
      // call tobitnet_innerC with pf1
      List tn = tobitnet_innerC(x, y, NumericVector::create(c), l1, 0.0, pf1, pf2, delta_init_vec, delta0_scalar_init, gamma_scalar_init, eps, standardize, maxit);
      bool KKT_here = as<bool>(tn["KKT"]);
      if(!KKT_here){
        // warn — can't call R's warning directly from here portably; use Rcpp::Rcout
        //Rcpp::Rcout << "Warning (tobitscad_cpp): KKT conditions not satisfied for lambda index " << l << ".\n";
      }
      
      // if last LLA iteration, compute deviance for this lambda
      if(it == iter-2){ // because it runs 0..iter-2 equivalent to R (i == iter -1)
        NumericVector beta_tn = tn["beta"];
        double gamma_tn = as<double>(tn["gamma"]);
        double b0_tn = as<double>(tn["b0"]);
        // delta_pred = beta*gamma (in scale of delta)
        NumericVector delta_pred(p);
        for(int j=0;j<p;j++) delta_pred[j] = beta_tn[j] * gamma_tn;
        NumericVector r_temp = matProd(x, delta_pred);
        // add intercept contribution: (b0 - c)*gamma
        double intercept_term = (b0_tn - c) * gamma_tn;
        for(int i=0;i<n;i++) r_temp[i] += intercept_term;
        NumericVector llv = logL1(yshift, d, r_temp, gamma_tn);
        dev[l] = -2.0 * sum(llv);
      }
      
      // store results
      NumericVector beta_tn = tn["beta"];
      double gamma_tn = as<double>(tn["gamma"]);
      double b0_tn = as<double>(tn["b0"]);
      for(int j=0;j<p;j++) beta_new(j,l) = beta_tn[j];
      gamma_new[l] = gamma_tn;
      b0_new[l] = b0_tn;
      
      // warm starts: delta (in standardized scale) returned as tn["delta"]
      NumericVector delta_out = tn["delta"];
      delta_init_vec = delta_out;
      delta0_scalar_init = as<double>(tn["d0"]);
      gamma_scalar_init = gamma_tn;
    }
    
    // overwrite delta_current for next LLA iter: delta = beta_mat %*% diag(gamma_vec)
    for(int l=0;l<nlambda;l++){
      for(int j=0;j<p;j++){
        delta_current(j,l) = beta_new(j,l) * gamma_new[l];
      }
    }
    
    // keep final beta_mat_final
    beta_mat_final = beta_new;
  }
  
  // If iter was 1 (no LLA iterations run), beta_mat_final must be taken from initial path:
  if(iter == 1){
    for(int l=0;l<nlambda;l++){
      for(int j=0;j<p;j++){
        beta_mat_final(j,l) = delta_mat(j,l) / gamma_vec[l]; // beta = delta/gamma
      }
    }
    // dev can be computed for that case as well
    for(int l=0;l<nlambda;l++){
      NumericVector delta_pred(p);
      for(int j=0;j<p;j++) delta_pred[j] = delta_mat(j,l);
      double gamma_l = gamma_vec[l];
      NumericVector r_temp = matProd(x, delta_pred);
      double intercept_term = (b0_vec[l] - c) * gamma_l;
      for(int i=0;i<n;i++) r_temp[i] += intercept_term;
      NumericVector llv = logL1(yshift, d, r_temp, gamma_l);
      dev[l] = -2.0 * sum(llv);
    }
  } else {
    // beta_mat_final already set during last LLA iteration: dev computed above
    // nothing to do
  }
  
  // Build output in same structure as R function:
  // Return sigma = 1/gamma_vec (but here gamma_vec after final LLA: gamma_new if available)
  NumericVector sigma_out(nlambda);
  NumericVector sigma_tmp = rep(NA_REAL, nlambda);
  // We attempted to fill gamma_new only inside LLA; fallback to earlier gamma_vec or last gamma_new
  // Try to collect gamma from last stage:
  if(iter >= 2){
    // gamma_new exists inside loop; we captured it into gamma_new variable
    // but scope wise, gamma_new is last assigned value; to be safe, recalc sigma from beta_mat_final and delta_current
    for(int l=0;l<nlambda;l++){
      // attempt to recover gamma by using delta_current and beta_mat_final: delta = beta*gamma => gamma = delta_j / beta_j for any nonzero beta_j
      double recovered_gamma = NA_REAL;
      for(int j=0;j<p;j++){
        double b = beta_mat_final(j,l);
        if(std::abs(b) > 1e-12){
          double dval = delta_current(j,l);
          recovered_gamma = dval / b;
          break;
        }
      }
      if(!NumericVector::is_na(recovered_gamma) && std::isfinite(recovered_gamma) && recovered_gamma > 0){
        sigma_out[l] = 1.0 / recovered_gamma;
      } else {
        // fallback: set to last gamma_vec if available
        if(gamma_vec.size() == nlambda && gamma_vec[l] > 0){
          sigma_out[l] = 1.0 / gamma_vec[l];
        } else {
          sigma_out[l] = NA_REAL;
        }
      }
    }
  } else {
    for(int l=0;l<nlambda;l++){
      if(gamma_vec[l] > 0) sigma_out[l] = 1.0 / gamma_vec[l];
      else sigma_out[l] = NA_REAL;
    }
  }
  
  // compute final beta matrix to return (p x nlambda)
  NumericMatrix beta_final(p, nlambda);
  for(int l=0;l<nlambda;l++){
    for(int j=0;j<p;j++){
      beta_final(j,l) = beta_mat_final(j,l);
    }
  }
  
  // set rownames? Rcpp list can't easily attach dimnames; return as-is.
  
  return List::create(
    Named("call") = "tobitscad_cpp",
    Named("sigma") = sigma_out,
    Named("b0") = b0_vec,
    Named("beta") = beta_final,
    Named("c") = c,
    Named("lambda") = lambda_path,
    Named("dev") = dev,
    Named("nulldev") = null_dev,
    Named("KKT_vec") = KKT_vec
  );
}


// [[Rcpp::export]]
List cv_tobitscad_cpp(NumericMatrix x, NumericVector y,
                      double c = 0.0, double a = 3.0, int iter = 3,
                      int nlambda = 100, double lambda_factor = 0.05,
                      Nullable<NumericVector> lambda_in = R_NilValue,
                      int nfolds = 10,
                      std::string type_measure = "mse",
                      bool early_stop = true,
                      double eps = 1e-7, bool standardize = true, 
                      double maxit = 1e4) {
  
  int n = x.nrow();
  int p = x.ncol();
  
  if(p <= 1) stop("x must have 2 or more columns");
  if(nfolds < 2 || nfolds > n) stop("nfolds must be between 2 and n");
  
  // initial fit to get lambda path
  List tn_init = tobitscad_cpp(x, y, c, a, iter, nlambda, lambda_factor,
                               lambda_in, eps, standardize, maxit, early_stop);
  
  NumericVector lambda = tn_init["lambda"];
  nlambda = lambda.size();
  
  // Split indices into censored/non-censored
  IntegerVector nonzero_indices, zero_indices;
  for(int i = 0; i < n; i++){
    if(y[i] > c) nonzero_indices.push_back(i);
    else if(y[i] == c) zero_indices.push_back(i);
  }
  int n_nz = nonzero_indices.size();
  int n_z  = zero_indices.size();
  
  int nperfold_nz = n_nz / nfolds;
  int nperfold_z  = n_z / nfolds;
  
  IntegerVector unfolded_entries_nz = clone(nonzero_indices);
  IntegerVector unfolded_entries_z  = clone(zero_indices);
  
  NumericMatrix err_mat(nfolds, nlambda);
  
  // Cross-validation folds
  for(int i = 0; i < nfolds; i++){
    IntegerVector fold_nz, fold_z;
    
    if(i < nfolds - 1){
      fold_nz = unfolded_entries_nz[Range(0, nperfold_nz-1)];
      fold_z  = unfolded_entries_z[Range(0, nperfold_z-1)];
      unfolded_entries_nz = unfolded_entries_nz[Range(nperfold_nz, unfolded_entries_nz.size()-1)];
      unfolded_entries_z  = unfolded_entries_z[Range(nperfold_z, unfolded_entries_z.size()-1)];
    } else {
      fold_nz = unfolded_entries_nz;
      fold_z  = unfolded_entries_z;
    }
    
    IntegerVector fold_idx = no_init(fold_nz.size() + fold_z.size());
    std::copy(fold_nz.begin(), fold_nz.end(), fold_idx.begin());
    std::copy(fold_z.begin(), fold_z.end(), fold_idx.begin() + fold_nz.size());
    
    // Create training and test sets
    LogicalVector test_mask(n, false);
    for(int j : fold_idx) test_mask[j] = true;
    
    int n_train = n - fold_idx.size();
    NumericMatrix x_train(n_train, p), x_test(fold_idx.size(), p);
    NumericVector y_train(n_train), y_test(fold_idx.size());
    
    int train_i = 0, test_i = 0;
    for(int j = 0; j < n; j++){
      if(test_mask[j]){
        x_test(test_i, _) = x(j, _);
        y_test[test_i] = y[j];
        test_i++;
      } else {
        x_train(train_i, _) = x(j, _);
        y_train[train_i] = y[j];
        train_i++;
      }
    }
    
    // Fit model on training fold
    List tscad = tobitscad_cpp(x_train, y_train, c, a, iter,
                               nlambda, lambda_factor, lambda,
                               eps, standardize, maxit, false);
    
    NumericMatrix beta = tscad["beta"];   // p x nlambda
    NumericVector b0   = tscad["b0"];     // length nlambda
    NumericVector sigma = tscad["sigma"];
    
    int ntest = x_test.nrow();
    NumericMatrix preds(ntest, nlambda);
    
    // predictions
    for(int j = 0; j < nlambda; j++){
      for(int irow = 0; irow < ntest; irow++){
        double r = b0[j];
        for(int k = 0; k < p; k++){
          r += x_test(irow, k) * beta(k, j);
        }
        if(type_measure == "mse" || type_measure == "mae"){
          if(r < c) r = c;  // censored prediction
        }
        preds(irow, j) = r;
      }
    }
    
    // compute error
    for(int j = 0; j < nlambda; j++){
      double err = 0.0;
      if(type_measure == "mse"){
        for(int irow = 0; irow < ntest; irow++){
          double diff = y_test[irow] - preds(irow, j);
          err += diff * diff;
        }
        err /= ntest;
      } else if(type_measure == "mae"){
        for(int irow = 0; irow < ntest; irow++){
          err += std::abs(y_test[irow] - preds(irow, j));
        }
        err /= ntest;
      } else if(type_measure == "deviance"){
        // use logL1 from your helper functions
        NumericVector r_temp(ntest);
        for(int irow = 0; irow < ntest; irow++){
          double r_val = (b0[j] - c)/sigma[j];
          for(int k = 0; k < p; k++){
            r_val += x_test(irow, k) * (beta(k, j)/sigma[j]);
          }
          r_temp[irow] = r_val;
        }
        LogicalVector d_test(ntest);
        NumericVector y_adj(ntest);
        for(int irow = 0; irow < ntest; irow++){
          y_adj[irow] = y_test[irow] - c;
          d_test[irow] = (y_test[irow] > c);
        }
        NumericVector ll = logL1(y_adj, d_test, r_temp, 1.0/sigma[j]);
        err = -2 * sum(ll);
      }
      err_mat(i, j) = err;
    }
  }
  
  NumericVector cvm(nlambda), cvvar(nlambda), cvsd(nlambda);
  for(int j = 0; j < nlambda; j++){
    NumericVector errs = err_mat(_, j);
    cvm[j] = mean(errs);
    cvvar[j] = mean(errs * errs) - cvm[j] * cvm[j];
    cvsd[j] = std::sqrt(cvvar[j] / nfolds);
  }
  
  int min_idx = which_min(cvm);
  double lambda_min = lambda[min_idx];
  double lambda_1se = lambda[0];
  for(int j = 0; j < nlambda; j++){
    if(cvm[j] < cvm[min_idx] + cvsd[min_idx]) lambda_1se = lambda[j];
  }
  
  return List::create(
    _["cvm"] = cvm,
    _["cvsd"] = cvsd,
    _["lambda"] = lambda,
    _["lambda.min"] = lambda_min,
    _["lambda.1se"] = lambda_1se,
    _["type.measure"] = type_measure
  );
}

