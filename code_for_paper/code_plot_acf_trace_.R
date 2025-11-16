# Required libraries
#devtools::install_github('TateJacobson/tobitnet')
library(tobitnet); 
library(truncnorm)  # for truncated normal sampling
library(MASS)       # for mvrnorm
library(Rcpp);sourceCpp('tobit_bayes.cpp'); sourceCpp('tobit_SnS.cpp')

# Simulate data
n <- 200
n_test <- 1000
n_all = n + n_test
p <- 300
s0 = 10
beta_true <- rep(0, p) ; beta_true[1:s0] <- c( rep(1, s0/2) , rep(-1, s0/2) )
rho <- 0.
Sigma <- outer(1:p, 1:p, function(i, j)rho^abs(i - j))
LL = chol(Sigma) 


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
  post_sams_hs <- result_cpp$beta_samples
  
  par(mar=c(2,2,1,1),mfrow=c(3,3))
  plot(post_sams_hs[,1],type = 'l',ylab = '',xlab = ''); abline(h=1,col='red',lwd=1)
  plot(post_sams_hs[,3],type = 'l',ylab = '',xlab = ''); abline(h=1,col='red',lwd=1)
  plot(post_sams_hs[,4],type = 'l',ylab = '',xlab = ''); abline(h=1,col='red',lwd=1)
  plot(post_sams_hs[,6],type = 'l',ylab = '',xlab = ''); abline(h=-1,col='red',lwd=1)
  plot(post_sams_hs[,7],type = 'l',ylab = '',xlab = ''); abline(h=-1,col='red',lwd=1)
  plot(post_sams_hs[,9],type = 'l',ylab = '',xlab = ''); abline(h=-1,col='red',lwd=1)
  plot(post_sams_hs[,11],type = 'l',ylab = '',xlab = ''); abline(h=0,col='red',lwd=1)
  plot(post_sams_hs[,15],type = 'l',ylab = '',xlab = ''); abline(h=0,col='red',lwd=1)
  plot(post_sams_hs[,21],type = 'l',ylab = '',xlab = ''); abline(h=0,col='red',lwd=1)
  
  library(coda)
  
  effectiveSize(post_sams_hs[,3])
  effectiveSize(post_sams_hs[,4])
  effectiveSize(post_sams_hs[,6])
  effectiveSize(post_sams_hs[,7])
  effectiveSize(post_sams_hs[,9])
  effectiveSize(post_sams_hs[,11])
  effectiveSize(post_sams_hs[,15])
  effectiveSize(post_sams_hs[,21])
  
  
  par(mar=c(2,2,2,1),mfrow=c(3,3))
  acf(post_sams_hs[,1],ylab = '',xlab = '',main= '' );
  title(main = paste( 'ESS=',round(effectiveSize(post_sams_hs[,1]),0)), line = 0.5) 
  acf(post_sams_hs[,3],main='',ylab = '',xlab = ''); 
  title(main = paste( 'ESS=',round(effectiveSize(post_sams_hs[,3]),0)), line = 0.5) 
  acf(post_sams_hs[,4],main='',ylab = '',xlab = '');
  title(main = paste( 'ESS=',round(effectiveSize(post_sams_hs[,4]),0)), line = 0.5) 
  acf(post_sams_hs[,6],main='',ylab = '',xlab = '');
  title(main = paste( 'ESS=',round(effectiveSize(post_sams_hs[,6]),0)), line = 0.5) 
  acf(post_sams_hs[,7],main='',ylab = '',xlab = '');
  title(main = paste( 'ESS=',round(effectiveSize(post_sams_hs[,7]),0)), line = 0.5) 
  acf(post_sams_hs[,9],main='',ylab = '',xlab = '');
  title(main = paste( 'ESS=',round(effectiveSize(post_sams_hs[,9]),0)), line = 0.5) 
  acf(post_sams_hs[,11],main='',ylab = '',xlab = '');
  title(main = paste( 'ESS=',round(effectiveSize(post_sams_hs[,11]),0)), line = 0.5) 
  acf(post_sams_hs[,15],main='',ylab = '',xlab = '');
  title(main = paste( 'ESS=',round(effectiveSize(post_sams_hs[,15]),0)), line = 0.5) 
  acf(post_sams_hs[,21],main='',ylab = '',xlab = '')
  title(main = paste( 'ESS=',round(effectiveSize(post_sams_hs[,21]),0)), line = 0.5) 
  
  par(mar=c(2,2,2,1),mfrow=c(1,2))
  plot(sam_hs_cpp$tau_sample[(10 + 1):5000, ] ,type = 'l',ylab = '',xlab = ''); abline(h=0,col='red',lwd=1)
  acf(sam_hs_cpp$tau_sample  ,main='',ylab = '',xlab = '')
  title(main = paste( 'ESS=',round(effectiveSize(sam_hs_cpp$tau_sample ),0)), line = 0.5) 
  
  


