# Required libraries
#devtools::install_github('TateJacobson/tobitnet')
library(tobitnet); library(tictoc)
library(truncnorm)  # for truncated normal sampling
library(MASS)       # for mvrnorm
library(Rcpp);sourceCpp('tobit_bayes.cpp'); sourceCpp('tobit_scad_fast.cpp')
sourceCpp('tobit_SnS.cpp')
#-------------------------------
# Tobit regression Gibbs sampler
#-------------------------------
source('Gibbs_tobit_hs.R')

# Simulate data
n <- 100
n_test <- 1000
n_all = n + n_test
p <- 120
s0 = 4
beta_true <- rep(0, p) ; beta_true[1:s0] <- c( rep(1, s0/2) , rep(-1, s0/2) )
rho <- 0.
Sigma <- outer(1:p, 1:p, function(i, j)rho^abs(i - j))
LL = chol(Sigma) 

HSout = tblassoout = scad_out = SnS_out = list()
for (ss in 1:100) {
  xall <- matrix(rnorm(n_all * p), n_all, p)%*% LL
  X <- xall[1:n,]
  xtest <- xall[-(1:n),]
  z_true <- xall %*% beta_true + rnorm(n_all ,sd = 1) # rnorm(n_all ,sd = 1)
  y = z_true[1:n] ; 
  censored <- 0 # quantile(y,probs = 1/4)
  y = pmax(censored , y ) 
  ytest = z_true[-(1:n)]
  ytest = pmax(censored , ytest ) 
  # Run Gibbs sampler
  
  result_cpp <- tobit_bayes_cpp(y, X, c_sensored = censored, n_iter = 5000, burn_in = 1000) 
  horseSH_cpp <- colMeans(result_cpp$beta_samples)
  
  res_SnS_cpp <- tobit_spikeslab_cpp(y, X, c_sensored = censored, n_iter = 5000, burn_in = 1000) 
  hat_SnS_cpp <- colMeans(res_SnS_cpp$beta_samples)
  
  tnet1_cv <- cv.tobitnet(x = X, y = y, c = censored, nfolds = 5)
  tb_lasso <- c(tobitnet(x = X, y = y, c = censored , lambda1 = tnet1_cv$lambda1.min)$beta) 
  tnet22_cv <- cv_tobitscad_cpp(x = X, y = y, c = censored, nlambda = 20, nfolds = 5)
  tb_scad2 <- c(tobitscad_cpp(x = X, y = y, c = censored , lambda_in = tnet22_cv$lambda.min)$beta) 
  
  ci_hs <- credible_intervals(result_cpp$beta_samples)
  ci_Sns <- credible_intervals(res_SnS_cpp$beta_samples)
  
  
  HSout[[ss]] <-  c(sum( (horseSH_cpp - beta_true)^2 ), mean((X%*%beta_true -X%*%horseSH_cpp )^2), mean((y -pmax(X%*%horseSH_cpp,censored) )^2),mean(( pmax(ytest,censored) -pmax(xtest%*%horseSH_cpp,censored) )^2),
                    selection_metrics(c(1:s0), selected_support = (1:p)[ !(ci_hs[, "Lower"] < 0 & ci_hs[, "Upper"] > 0)] , p = p))
  scad_out[[ss]] <-   c(sum( (tb_scad2 - beta_true)^2 ), mean((X%*%beta_true -X%*%tb_scad2 )^2), mean((y -pmax(X%*%tb_scad2,censored) )^2),mean((pmax(ytest,censored)  -pmax(xtest%*%tb_scad2,censored) )^2),
                        selection_metrics(true_support = 1:s0,selected_support = which(tb_scad2 !=0),p = p) )
  tblassoout[[ss]] <- c(sum( (tb_lasso - beta_true)^2 ), mean((X%*%beta_true -X%*%tb_lasso )^2), mean((y -pmax(X%*%tb_lasso,0) )^2),mean((ytest -pmax(xtest%*%tb_lasso,censored) )^2),
                        selection_metrics(true_support = 1:s0,selected_support = which(tb_lasso !=0),p = p)  )
  SnS_out[[ss]] <-  c(sum( (hat_SnS_cpp - beta_true)^2 ), mean((X%*%beta_true -X%*%hat_SnS_cpp )^2), mean((y -pmax(X%*%hat_SnS_cpp,censored) )^2),mean(( pmax(ytest,censored) -pmax(xtest%*%hat_SnS_cpp,censored) )^2),
                      selection_metrics(c(1:s0), selected_support = (1:p)[ !(ci_Sns[, "Lower"] < 0 & ci_Sns[, "Upper"] > 0)] , p = p) )
  
  print(ss)
}
rm(Sigma,LL,X,xall)
save.image(file = 'tobit_p120n100_s04_.rda')






