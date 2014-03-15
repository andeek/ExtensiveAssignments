### Code for MCMC for Regional Model ###
#source('/code/library.R')
#source('/code/data_format.R')
#source('/code/helpers.R')
#load('/data/res1.rda)
#load('/data/res1.rda)
#load('/data/res1.rda)

#### Conditional Distributions ####
cond_b0 <- function(beta0, beta1, beta2, beta3, sigma0, x, y){
  delta <- y==0
  p <- exp(beta0 + beta1*x)/(1+exp(beta0 + beta1*x))
  lambda <- exp(beta2 + beta3*x)
  sum(delta*log(1-p+p*exp(-lambda)) + (1-delta)*log(p)) - 1/2*(beta0)^2/sigma0^2
}

cond_b1 <- function(beta0, beta1, beta2, beta3, sigma1, x, y){
  delta <- y==0
  p<-exp(beta0 + beta1*x)/(1+exp(beta0 + beta1*x))
  lambda<-exp(beta2 + beta3*x)
  sum(delta*log(1-p+p*exp(-lambda)) + (1-delta)*log(p)) - 1/2*(beta1)^2/sigma1^2
}

cond_b2 <- function(beta0, beta1, beta2, beta3, sigma2, x, y){
  delta <- y==0
  p<-exp(beta0 + beta1*x)/(1+exp(beta0 + beta1*x))
  lambda<-exp(beta2 + beta3*x)
  sum(delta*log(1-p+p*exp(-lambda)) + (1-delta)*(y*(beta2+beta3*x) - lambda)) - 1/2*(beta2)^2/sigma2^2
}

cond_b3 <- function(beta0, beta1, beta2, beta3, sigma3, x, y){
  delta <- y==0
  p<-exp(beta0 + beta1*x)/(1+exp(beta0 + beta1*x))
  lambda<-exp(beta2 + beta3*x)
  sum(delta*log(1-p+p*exp(-lambda)) + (1-delta)*(y*(beta2+beta3*x) - lambda)) - 1/2*(beta3)^2/sigma3^2
}

#### Sampling Functions ####
sample_beta0 <- function(y, x, beta0, beta1, beta2, beta3, sigma0, sigma00, b0) {
  beta0_star <- rnorm(1, beta0, sigma00)
  
  if (beta0_star < -b0 | beta0_star > b0) return(beta0)
  
  lrho = cond_b0(beta0=beta0_star, beta1, beta2, beta3, sigma0, x, y) - 
    cond_b0(beta0=beta0, beta1, beta2, beta3, sigma0, x, y)
  
  return(ifelse(log(runif(1))<lrho, beta0_star, beta0))
}

sample_beta1 <- function(y, x, beta0, beta1, beta2, beta3, sigma1, sigma01, b1) {
  beta1_star <- rnorm(1, beta1, sigma01)
  
  if (beta1_star < -b1 | beta1_star > b1) return(beta1)
  
  lrho <- cond_b1(beta0, beta1=beta1_star, beta2, beta3, sigma1, x, y) - 
    cond_b1(beta0, beta1=beta1, beta2, beta3, sigma1, x, y)
  
  return(ifelse(log(runif(1))<lrho, beta1_star, beta1))
}

sample_beta2 <- function(y, x, beta0, beta1, beta2, beta3, sigma2, sigma02, b2) {
  beta2_star <- rnorm(1, beta2, sigma02)
  
  if (beta2_star < -b2 | beta2_star > b2) return(beta2)
  
  lrho <- cond_b2(beta0, beta1, beta2=beta2_star, beta3, sigma2, x, y) - 
    cond_b2(beta0, beta1, beta2=beta2, beta3, sigma2, x, y)
  
  return(ifelse(log(runif(1))<lrho, beta2_star, beta2))
}

sample_beta3 <- function(y, x, beta0, beta1, beta2, beta3, sigma3, sigma03, b3) {
  beta3_star <- rnorm(1, beta3, sigma03)
  
  if (beta3_star < -b3 | beta3_star > b3) return(beta3)
  
  lrho <- cond_b3(beta0, beta1, beta2, beta3=beta3_star, sigma3, x, y) - 
    cond_b3(beta0, beta1, beta2, beta3=beta3, sigma3, x, y)
  
  return(ifelse(log(runif(1))<lrho, beta3_star, beta3))
}

sample_sigma<-function(z, beta, a0, g0){
  a<-length(z)*(a0+3/2)-1
  b<-sum(beta^2/2) + g0
  sqrt(1/rgamma(1, a, b))
}

run_mcmc = function(dat, beta0, beta1, beta2, beta3,
                            sigma0, sigma1, sigma2, sigma3, s=1, 
                            n.reps=10, b0=100, b1=100, b2=100, b3=100, a=.5, b=.5, 
                            counter=TRUE, tune=TRUE){
  z=unique(dat$store)
  sigma00 = sigma01 = sigma02 = sigma03 = s
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
  sigmat_keep <- matrix(numeric(n.reps), nrow=n.reps, ncol=4)
  sigma_keep[1,] <- c(sigma00, sigma01, sigma02, sigma03)
  
  for (i in 2:n.reps) {      
    for(j in 1:length(z)){
      y<-dat[dat$store==z[j],]$mvm
      x<-dat[dat$store==z[j],]$price
      beta0_old<-beta0[j]
      beta1_old<-beta1[j]
      beta2_old<-beta2[j]
      beta3_old<-beta3[j]
      
      beta0[j] <- sample_beta0(y, x, beta0[j], beta1[j], beta2[j], beta3[j], sigma0, sigma00, b0)
      beta1[j] <- sample_beta1(y, x, beta0[j], beta1[j], beta2[j], beta3[j], sigma1, sigma01, b1)
      beta2[j] <- sample_beta2(y, x, beta0[j], beta1[j], beta2[j], beta3[j], sigma2, sigma02, b2)
      beta3[j] <- sample_beta3(y, x, beta0[j], beta1[j], beta2[j], beta3[j], sigma3, sigma03, b3)
    
      if (tune) {
        if(beta0[j]==beta0_old){
          sigma00 = sigma00/1.1 
        }else{
          sigma00 = sigma00*1.1
        }
        
        if(beta1[j]==beta1_old){
          sigma01 = sigma01/1.1 
        }else{
          sigma01 = sigma01*1.1
        }
        
        if(beta2[j]==beta2_old){
          sigma02 = sigma02/1.1 
        }else{
          sigma02 = sigma02*1.1
        }
        
        if(beta3[j]==beta3_old){
          sigma03 = sigma03/1.1 
        }else{
          sigma03 = sigma03*1.1
        }
      }
    }
    
    sigma0 <- sample_sigma(z, beta0, a, b)
    sigma1 <- sample_sigma(z, beta1, a, b)
    sigma2 <- sample_sigma(z, beta2, a, b)
    sigma3 <- sample_sigma(z, beta3, a, b)

    beta0_keep[i,] = beta0
    beta1_keep[i,] = beta1
    beta2_keep[i,] = beta2
    beta3_keep[i,] = beta3
    sigma_keep[i,] = c(sigma0, sigma1, sigma2, sigma3)
    sigmat_keep[i,] = c(sigma00, sigma01, sigma02, sigma03)
    
    if(counter){
      cat("\r")
      cat("Iter:", i, "\r")
    }
  }
  return(list(beta0=beta0_keep, beta1=beta1_keep, beta2=beta2_keep, beta3=beta3_keep, sigma=sigma_keep, sigmat=sigmat_keep))
}



#### Run 3 Chains ####
res1 <- run_mcmc(dat=dat, beta0=runif(length(unique(dat$store)),5,8), beta1=runif(length(unique(dat$store)),-12,-8), beta2=runif(length(unique(dat$store)),5,8), beta3=runif(length(unique(dat$store)),-10,-8),
                sigma0=1, sigma1=1, sigma2=1, sigma3=1, s=1,
                n.reps=3000, a=0.5, b=0.5)

res2 <- run_mcmc(dat=dat, beta0=runif(length(unique(dat$store)),5,8), beta1=runif(length(unique(dat$store)),-12,-8), beta2=runif(length(unique(dat$store)),5,8), beta3=runif(length(unique(dat$store)),-10,-8),
                sigma0=1, sigma1=6, sigma2=3, sigma3=2, s=1,
                n.reps=3000, a=0.5, b=0.5)

res3 <- run_mcmc(dat=dat, beta0=runif(length(unique(dat$store)),5,8), beta1=runif(length(unique(dat$store)),-12,-8), beta2=runif(length(unique(dat$store)),5,8), beta3=runif(length(unique(dat$store)),-10,-8),
                sigma0=6, sigma1=5, sigma2=8, sigma3=10, s=1,
                n.reps=3000, a=0.5, b=0.5)

#### Get Estimates and CI's ####
#### Median and CI's for each store ####
beta0<-rbind(res1$beta0[1000:3000,], res2$beta0[1000:3000,], res3$beta0[1000:3000,])
meds_b0<-apply(beta0, 2, median)
cis_b0<-apply(beta0, 2, function(u) sort(u)[c(150,5850)])

beta1<-rbind(res1$beta1[1000:3000,], res2$beta1[1000:3000,], res3$beta1[1000:3000,])
meds_b1<-apply(beta1, 2, median)
cis_b1<-apply(beta1, 2, function(u) sort(u)[c(150,5850)])

beta2<-rbind(res1$beta2[1000:3000,], res2$beta2[1000:3000,], res3$beta2[1000:3000,])
meds_b2<-apply(beta2, 2, median)
cis_b2<-apply(beta2, 2, function(u) sort(u)[c(150,5850)])

beta3<-rbind(res1$beta3[1000:3000,], res2$beta3[1000:3000,], res3$beta3[1000:3000,])
meds_b3<-apply(beta3, 2, median)
cis_b3<-apply(beta3, 2, function(u) sort(u)[c(150,5850)])

est.bayes<-cbind(b0=meds_b0, b1=meds_b1, b2=meds_b2, b3=meds_b3)
row.names(est.bayes)<-unique(dat$store)
cis.bayes<-cbind(t(cis_b0), t(cis_b1), t(cis_b2), t(cis_b3))
row.names(cis.bayes)<-unique(dat$store)
names(cis.bayes) <- c("b0_l", "b0_u","b1_l", "b1_u","b2_l", "b2_u","b3_l", "b3_u")
