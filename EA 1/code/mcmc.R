#### Libraries ####
library(coda)
library(plyr)
library(lubridate)
library(ggplot2)
library(reshape)

setwd("./GitHub/ExtensiveAssignments")

#### Conditional Distributions ####
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

#### Sampling Functions ####
sample_beta0 <- function(y, x, beta0, beta1, beta2, beta3, sigma0, sigma, b0) {
  beta0_star <- rnorm(1, beta0, sigma)
  
  if (beta0_star < -b0 | beta0_star > b0) return(beta0)
  
  lrho = cond_b0(beta0=beta0_star, beta1, beta2, beta3, sigma0, x, y) - 
    cond_b0(beta0=beta0, beta1, beta2, beta3, sigma0, x, y)
  
  return(ifelse(log(runif(1))<lrho, beta0_star, beta0))
}

sample_beta1 <- function(y, x, beta0, beta1, beta2, beta3, sigma1, sigma, b1) {
  beta1_star <- rnorm(1, beta1, sigma)
  
  if (beta1_star < -b1 | beta1_star > b1) return(beta1)
  
  lrho <- cond_b1(beta0, beta1=beta1_star, beta2, beta3, sigma1, x, y) - 
    cond_b1(beta0, beta1=beta1, beta2, beta3, sigma1, x, y)
  
  return(ifelse(log(runif(1))<lrho, beta1_star, beta1))
}

sample_beta2 <- function(y, x, beta0, beta1, beta2, beta3, sigma2, sigma, b2) {
  beta2_star <- rnorm(1, beta2, sigma)
  
  if (beta2_star < -b2 | beta2_star > b2) return(beta2)
  
  lrho <- cond_b2(beta0, beta1, beta2=beta2_star, beta3, sigma2, x, y) - 
    cond_b2(beta0, beta1, beta2=beta2, beta3, sigma2, x, y)
  
  return(ifelse(log(runif(1))<lrho, beta2_star, beta2))
}

sample_beta3 <- function(y, x, beta0, beta1, beta2, beta3, sigma3, sigma, b3) {
  beta3_star <- rnorm(1, beta3, sigma)
  
  if (beta3_star < -b3 | beta3_star > b3) return(beta3)
  
  lrho <- cond_b3(beta0, beta1, beta2, beta3=beta3_star, sigma3, x, y) - 
    cond_b3(beta0, beta1, beta2, beta3=beta3, sigma3, x, y)
  
  return(ifelse(log(runif(1))<lrho, beta3_star, beta3))
}

sample_sigma<-function(dat, beta, a0, g0){
  z<-unique(dat$store)
  y_stand<-c()
  for(i in 1:length(z)){
    y_stand<-c(y_stand, (dat[dat$store==z[i],]$mvm - beta[i])^2)
  }
  a<-length(dat$mvm)*(a0+1)+1
  b<-sum(y_stand)/2 + g0
  
  sqrt(1/rgamma(1, a, b))
}


run_mcmc = function(dat, beta0, beta1, beta2, beta3,
                            sigma0, sigma1, sigma2, sigma3, sigma=100,
                            n.reps=10, b0=100, b1=100, b2=100, b3=100, a=.5, b=.5, counter=TRUE){
  z=unique(dat$store)
  beta0_keep <- matrix(numeric(n.reps), nrow=n.reps, ncol=length(z))
  beta0_keep[1,] <- beta0
  beta1_keep <- matrix(numeric(n.reps), nrow=n.reps, ncol=length(z))
  beta1_keep[1,] <- beta1
  beta2_keep <- matrix(numeric(n.reps), nrow=n.reps, ncol=length(z))
  beta2_keep[1,] <- beta2
  beta3_keep <- matrix(numeric(n.reps), nrow=n.reps, ncol=length(z))
  beta3_keep[1,] <- beta3
  sigma_keep <- matrix(numeric(n.reps), nrow=n.reps, ncol=4)
  sigma_keep[1,] <- c(sigma0, sigma1, sigma2, sigma3)
  
  for (i in 1:n.reps) {      
    for(j in 1:length(z)){
      y<-dat[dat$store==z[j],]$mvm
      x<-dat[dat$store==z[j],]$price
      
      beta0[j] <- sample_beta0(y, x, beta0[j], beta1[j], beta2[j], beta3[j], sigma0, sigma, b0)
      beta1[j] <- sample_beta1(y, x, beta0[j], beta1[j], beta2[j], beta3[j], sigma1, sigma, b1)
      beta2[j] <- sample_beta2(y, x, beta0[j], beta1[j], beta2[j], beta3[j], sigma2, sigma, b2)
      beta3[j] <- sample_beta3(y, x, beta0[j], beta1[j], beta2[j], beta3[j], sigma3, sigma, b3)
    }
    
    sigma0 <- sample_sigma(dat, beta0, a, b)
    sigma1 <- sample_sigma(dat, beta1, a, b)
    sigma2 <- sample_sigma(dat, beta2, a, b)
    sigma3 <- sample_sigma(dat, beta3, a, b)

    
    beta0_keep[i,] = beta0
    beta1_keep[i,] = beta1
    beta2_keep[i,] = beta2
    beta3_keep[i,] = beta3
    sigma_keep[i,] = c(sigma0, sigma1, sigma2, sigma3)
    if(counter){
      cat("\r")
      cat("Iter:", i, "\r")
    }
  }
  return(list(beta0=beta0_keep, beta1=beta1_keep, beta2=beta2_keep, beta3=beta3_keep, sigma=sigma_keep))
}




#### Data ####
dat<-read.table(file='EA 1/data/greenbeandat.txt', header=TRUE)[,c(1,3,4,5)]
dat$date <-  ymd(as.character(dat$date))

#### Sample Stores ####
samp_store<-function(n, seed, dat){
  set.seed(seed)
  s<-sample(unique(dat$store),n)
  dat[dat$store %in% s,]
}

#### Plot Price vs. Volume ####
## Sample 4 stores
stores4<-samp_store(4, 2, dat=dat)

set.seed(1)
ptm <- proc.time()
res1 = run_mcmc(dat=stores4, beta0=runif(4,-10,10), beta1=runif(4,-10,10), beta2=runif(4,-10,10), beta3=runif(4,-10,10),
                sigma0=10, sigma1=10, sigma2=10, sigma3=10, sigma=0.5,
                        n.reps=500, a=0.5, b=0.5)
proc.time() - ptm

#### Check some acceptance probabilities? ####
sum(diff(res1$beta3[,1])!=0)/500
sum(diff(res1$beta2[,1])!=0)/500
sum(diff(res1$beta1[,1])!=0)/500
sum(diff(res1$beta0[,4])!=0)/500

qplot(1:500, res1$beta0[,1], geom="line", xlab="iteration", ylab=expression(beta[paste("0,1")]))
qplot(1:500, res1$beta1[,1], geom="line", xlab="iteration", ylab=expression(beta[paste("0,1")]))
qplot(1:500, res1$beta2[,1], geom="line", xlab="iteration", ylab=expression(beta[paste("0,1")]))
qplot(1:500, res1$beta3[,1], geom="line", xlab="iteration", ylab=expression(beta[paste("0,1")]))





res2 = run_mcmc(alpha0=runif(1,0,20), beta0=runif(1,0,2))
res3 = run_mcmc(alpha0=runif(1,0,20), beta0=runif(1,0,2))

#### OLD ####
run_mcmc_nosigma = function(dat, beta0, beta1, beta2, beta3, s=10, n.reps=10, tune=TRUE, b0=100, b1=100, b2=100, b3=100) {
  z=unique(dat$store)
  sigma0 = sigma1 = sigma2 = sigma3 = s
  beta0_keep <- matrix(numeric(n.reps), nrow=n.reps, ncol=length(z))
  beta0_keep[1,] <- beta0
  beta1_keep <- matrix(numeric(n.reps), nrow=n.reps, ncol=length(z))
  beta1_keep[1,] <- beta1
  beta2_keep <- matrix(numeric(n.reps), nrow=n.reps, ncol=length(z))
  beta2_keep[1,] <- beta2
  beta3_keep <- matrix(numeric(n.reps), nrow=n.reps, ncol=length(z))
  beta3_keep[1,] <- beta3
  
  for (i in 1:n.reps) {
    
    # Automatically tune alpha
    beta0_old = beta0
    beta1_old = beta1
    beta2_old = beta2
    beta3_old = beta3
    
    for(j in 1:length(unique(z))){
      y<-dat[dat$store==z[j],]$y
      x<-dat[dat$store==z[j],]$x
      
      beta0[j] <- sample_beta0(y, x, beta0[j], beta1[j], beta2[j], beta3[j], sigma0, b0)
      beta1[j] <- sample_beta1(y, x, beta0[j], beta1[j], beta2[j], beta3[j], sigma1, b1)
      beta2[j] <- sample_beta2(y, x, beta0[j], beta1[j], beta2[j], beta3[j], sigma2, b2)
      beta3[j] <- sample_beta3(y, x, beta0[j], beta1[j], beta2[j], beta3[j], sigma3, b3)
    }
    
    if (tune) {
      if (all(beta0==beta0_old)) {
        sigma0 = sigma0/1.1 
      } else {
        sigma0 = sigma0*1.1
      }
      if (all(beta1==beta1_old)) {
        sigma1 = sigma1/1.1 
      } else {
        sigma1 = sigma1*1.1
      }
      if (all(beta2==beta2_old)) {
        sigma2 = sigma2/1.1 
      } else {
        sigma2 = sigma2*1.1
      }
      if (all(beta3==beta3_old)) {
        sigma3 = sigma3/1.1 
      } else {
        sigma3 = sigma3*1.1
      }
    }
    
    beta0_keep[i,] = beta0
    beta1_keep[i,] = beta1
    beta2_keep[i,] = beta2
    beta3_keep[i,] = beta3
  }
  return(list(beta0=beta0_keep, beta1=beta1_keep, beta2=beta2_keep, beta3=beta3_keep, sigma=c(sigma0, sigma1, sigma2, sigma3)))
}




