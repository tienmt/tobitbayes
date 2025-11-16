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




tobit_bayes_WI <- function(y, X, c_sensored = 0,
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
    
    # Step 2: Sample beta (Woodbury version, stable)
    D <- diag(tau2 * lambda2)
    XD <- X %*% D
    
    # Compute S = sigma2 * I_n + X D X^T
    S <- sigma2 * diag(n) + XD %*% tX
    
    # Numerically stable inversion via Cholesky
    S_inv <- chol2inv(chol(S))
    
    # Posterior mean
    mu_beta <- D %*% tX %*% (S_inv %*% z)
    
    # Matrix for reuse
    B <- D %*% tX %*% S_inv
    
    # Random draws
    u <- rnorm(p) * sqrt(diag(D))       # u ~ N(0, D)
    v <- rnorm(n, sd = sqrt(sigma2))    # v ~ N(0, sigma^2 I)
    
    # Sample beta
    beta <- mu_beta + (diag(p) - B %*% X) %*% u + B %*% v
    
    
    
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







tobit_bayes_backCHol <- function(y, X, c_sensored = 0,
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
    
    # Step 2: Sample beta (Standard p x p version, fast if p < n)
    
    # 1. Calculate posterior precision
    # D is diag(tau2 * lambda2), so D_inv is diag(1 / (tau2 * lambda2))
    D_inv <- diag(1 / (tau2 * lambda2), nrow = p, ncol = p)
    post_precision <- (1/sigma2) * tXX + D_inv
    
    # 2. Calculate posterior mean component
    # We want mu = (post_precision^-1) * (1/sigma2 * tX %*% z)
    b <- (1/sigma2) * tX %*% z
    
    # 3. Use Cholesky for stable inversion and sampling
    # L is the *upper* triangle Cholesky factor: L^T L = post_precision
    L_post <- chol(post_precision) 
    
    # 4. Solve for mean: L^T L mu = b
    # This is more stable than mu = solve(post_precision, b)
    mu_beta <- backsolve(L_post, forwardsolve(t(L_post), b))
    
    # 5. Sample from N(0, post_precision^-1)
    # We solve L^T beta_dev = w, where w ~ N(0, I_p)
    # The resulting beta_dev ~ N(0, (L^T L)^-1) = N(0, post_precision^-1)
    w <- rnorm(p)
    beta_dev <- backsolve(L_post, w)
    
    # 6. Combine mean and deviation for the final sample
    beta <- mu_beta + beta_dev
    
    
    
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


# Simulate data
n <- 100
n_test <- 1000
n_all = n + n_test
p <- 300
xall <- matrix(rnorm(n_all * p), n_all, p)
X <- xall[1:n,]
xtest <- xall[-(1:n),]
s0 = 4
beta_true <- rep(0, p) ; 
beta_true[1:s0] <- c( rep(1, s0/2) , rep(-1, s0/2) )
z_true <- xall %*% beta_true + rnorm(n_all ,sd = 1)

y = z_true[1:n] ; 
censored <- quantile(y,probs = 1/3)
y = pmax(censored , y ) 

library(tictoc)

# fit the Gibbs sampler
tic(); res <- tobit_bayes(y, X,n_iter = 5000,  burn_in = 1000); toc()
tic(); res2 <- tobit_bayes_WI(y, X,n_iter = 5000,  burn_in = 1000); toc()
tic(); res3 <- tobit_bayes_backCHol(y, X,  n_iter = 5000,  burn_in = 1000); toc()

colMeans(res$beta_samples) [1:s0]
colMeans(res2$beta_samples) [1:s0]
colMeans(res3$beta_samples) [1:s0]

