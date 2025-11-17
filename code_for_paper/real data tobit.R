

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
library(abess)
data("trim32"); rownames(trim32) <- NULL
n = nrow(trim32)
p = 500
trim32$y <- scale(trim32$y)
trim32[,-1] <- scale(trim32[,-1])

censored = 0
HSout = tblassoout = scad_out = SnS_out = list()

for (ss in 1:100) {
  test = sample(1:n, 36)
  Y = trim32$y[-test]
  y = pmax(Y,0)
  X = as.matrix(trim32[-test,-1]) ; tX = t(X)
  Ytest =  trim32$y[test]
  Xtest = as.matrix(trim32[test,-1] )
  
  result_cpp <- tobit_bayes_cpp(y, X, c_sensored = censored, n_iter = 5000, burn_in = 1000) 
  horseSH_cpp <- colMeans(result_cpp$beta_samples)
  
  res_SnS_cpp <- tobit_spikeslab_cpp(y, X, c_sensored = censored, n_iter = 5000, burn_in = 1000) 
  hat_SnS_cpp <- colMeans(res_SnS_cpp$beta_samples)
  
  tnet1_cv <- cv.tobitnet(x = X, y = y, c = censored, nfolds = 5)
  tb_lasso <- c(tobitnet(x = X, y = y, c = censored , lambda1 = tnet1_cv$lambda1.min)$beta) 
  tnet22_cv <- cv_tobitscad_cpp(x = X, y = y, c = censored, nlambda = 20, nfolds = 5)
  tb_scad2 <- c(tobitscad_cpp(x = X, y = y, c = censored , lambda_in = tnet22_cv$lambda.min)$beta) 
  
  
  
  tblassoout[[ss]] = mean((pmax(Ytest,0)  -pmax(Xtest%*%tb_lasso,0) )^2)
  scad_out[[ss]] = mean((pmax(Ytest,0)  -pmax(Xtest%*%tb_scad2,0) )^2)
  HSout[[ss]] = mean((pmax(Ytest,0)  -pmax(Xtest%*%horseSH_cpp,0) )^2)
  SnS_out[[ss]] = mean((pmax(Ytest,0)  -pmax(Xtest%*%hat_SnS_cpp,0) )^2)
  print(ss)
}
save.image('out_realdata.rda')
mean(sapply(tblassoout, function(x) x[1])); sd(sapply(tblassoout, function(x) x[1]))
mean(sapply(scad_out, function(x) x[1])); sd(sapply(scad_out, function(x) x[1]))
mean(sapply(HSout, function(x) x[1])); sd(sapply(HSout, function(x) x[1]))
mean(sapply(SnS_out, function(x) x[1])); sd(sapply(SnS_out, function(x) x[1]))


