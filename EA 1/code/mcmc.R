library(coda)
library(plyr)

cond_b0 <- function(beta0, beta1, beta2, beta3, sigma0, x, y){
  delta <- y==0
  p<-exp(beta0 + beta1*x)/(1+exp(beta0 + beta1*x))
  lambda<-exp(beta2 + beta3*x)
  sum(delta*log(1-p+p*exp(-lambda)) + (1-delta)*log(p) - 1/2*(y-beta0)^2/sigma0^2)
}

cond_b1 <- function(beta0, beta1, beta2, beta3, sigma1, x, y){
  delta <- y==0
  p<-exp(beta0 + beta1*x)/(1+exp(beta0 + beta1*x))
  lambda<-exp(beta2 + beta3*x)
  sum(delta*log(1-p+p*exp(-lambda)) + (1-delta)*log(p) - 1/2*(y-beta1)^2/sigma1^2)
}

cond_b2 <- function(beta0, beta1, beta2, beta3, sigma2, x, y){
  delta <- y==0
  p<-exp(beta0 + beta1*x)/(1+exp(beta0 + beta1*x))
  lambda<-exp(beta2 + beta3*x)
  sum(delta*log(1-p+p*exp(-lambda)) + (1-delta)*(y*(beta2+beta3*x) - lambda) - 1/2*(y-beta2)^2/sigma2^2)
}

cond_b3 <- function(beta0, beta1, beta2, beta3, sigma3, x, y){
  delta <- y==0
  p<-exp(beta0 + beta1*x)/(1+exp(beta0 + beta1*x))
  lambda<-exp(beta2 + beta3*x)
  sum(delta*log(1-p+p*exp(-lambda)) + (1-delta)*(y*(beta2+beta3*x) - lambda) - 1/2*(y-beta3)^2/sigma3^2)
}

cond_sigma <- function(beta, sigma, x, y, nu, tau){
  sum(-1/2*(y-beta)^2/sigma^2 - nu*tau^2/(2*sigma^2) - (1+nu/2)*log(sigma^2))
}

sample_beta0 <- function(y, x, beta0, beta1, beta2, beta3, sigma0, b0) {
  beta0_star <- rnorm(1, beta0, sigma0)
  
  if (beta0_star < 0 | beta0_star > b0) return(beta0)
  
  lrho = cond_bo(beta0=beta0_star, beta1, beta2, beta3, sigma0, x, y) - 
    cond_bo(beta0=beta0, beta1, beta2, beta3, sigma0, x, y)
  
  return(ifelse(log(runif(1))<lrho, beta0_star, beta0))
}

sample_beta0 <- function(y, x, beta0, beta1, beta2, beta3, sigma1, b1) {
  beta1_star <- rnorm(1, beta1, sigma1)
  
  if (beta0_star < 0 | beta0_star > b1) return(beta1)
  
  lrho <- cond_b1(beta0, beta1=beta1_star, beta2, beta3, sigma1, x, y) - 
    cond_bo(beta0, beta1=beta1, beta2, beta3, sigma1, x, y)
  
  return(ifelse(log(runif(1))<lrho, beta1_star, beta1))
}

sample_beta2 <- function(y, x, beta0, beta1, beta2, beta3, sigma2, b2) {
  beta2_star <- rnorm(1, beta2, sigma2)
  
  if (beta2_star < 0 | beta2_star > b2) return(beta2)
  
  lrho <- cond_bo(beta0, beta1, beta2=beta2_star, beta3, sigma2, x, y) - 
    cond_bo(beta0, beta1, beta2=beta2, beta3, sigma0, x, y)
  
  return(ifelse(log(runif(1))<lrho, beta2_star, beta2))
}


sample_beta3 <- function(y, x, beta0, beta1, beta2, beta3, sigma3, b3) {
  beta3_star <- rnorm(1, beta3, sigma3)
  
  if (beta3_star < 0 | beta3_star > b3) return(beta3)
  
  lrho <- cond_bo(beta0, beta1, beta2, beta3=beta3_star, sigma0, x, y) - 
    cond_bo(beta0, beta1, beta2, beta3=beta3, sigma0, x, y)
  
  return(ifelse(log(runif(1))<lrho, beta3_star, beta3))
}

sample_sigma<-function(y, x, beta, nu, tau, sigma, o1)
  sigma_star <- 1/rchisq(1, nu)

  if (sigma_star < 0 | sigma_star > o1) return(sigma)

  lrho <- cond_bo(beta, sigma=sigma_star, x, y, nu, tau) - 
    cond_bo(beta, sigma=sigma, x, y, nu, tau)

  return(ifelse(log(runif(1))<lrho, beta3_star, beta3))
}
