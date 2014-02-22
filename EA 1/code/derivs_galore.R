#### Libraries ####
library(lubridate)

#### Function for newtraph ####
ders.bean <- function(ps, dat) {
  ##likelihood
  b0 <- ps[1]
  b1 <- ps[2]
  b2 <- ps[3]
  b3 <- ps[4]  
  
  x <- dat$price
  y <- dat$mvm
  
  T1 <- exp(b0 + b1*x)
  T2 <- exp(b2 + b3*x)  
  
  p <- T1/(1+T1)
  lambda <- T2
  
  delta<-rep(0,length(y))
  delta[y==0]<-1

  loglik <- sum(delta*log(1-p + p*exp(-lambda)) + (1-delta)*(log(p) - lfactorial(y) + y*log(lambda) - lambda))
  
  ##gradient
  dldp <- delta*(exp(-lambda)-1)/(1-p + p*exp(-lambda)) + (1-delta)/p 
  dldlam <- delta*(-p*exp(-lambda))/(1-p + p*exp(-lambda)) + (1-delta)*(y/lambda-1)
  dpdb0 <- T1/(1+T1)^2
  dpdb1 <- x*T1/(1+T1)^2
  dpdb2 <- 0
  dpdb3 <- 0
  dlamdb0 <- 0
  dlamdb1 <- 0
  dlamdb2 <- T2
  dlamdb3 <- x*T2
   
  grad <- c(sum(dldp*dpdb0), sum(dldp*dpdb1), sum(dldlam*dlamdb2), sum(dldlam*dlamdb3))
   
  ##hessian
  d2ldp2 <- -delta*((exp(-lambda)-1)/(1-p + p*exp(-lambda)))^2 - (1-delta)/p^2 
  d2ldpdlam <- -delta*(exp(-lambda))/((1-p) + p*exp(-lambda))^2
  d2ldlam2 <- delta*p*(1-p*exp(-lambda))/(1-p + p*exp(-lambda))^2 - (1-delta)*(y/lambda^2)
  d2pdb02 <- (T1/(1+T1))^2 - 2*(T1/(1+T1))^3
  d2pdb12 <- x^2*((T1/(1+T1))^2 - 2*(T1/(1+T1))^3)
  d2pdb22 <- 0
  d2pdb32 <- 0  
  d2pdb0db1 <- x*((T1/(1+T1))^2 - 2*(T1/(1+T1))^3)
  d2pdb0db2 <- 0
  d2pdb0db3 <- 0
  d2pdb1db2 <- 0
  d2pdb1db3 <- 0
  d2pdb2db3 <- 0

  d2lamdb02 <- 0
  d2lamdb12 <- 0
  d2lamdb22 <- T2
  d2lamdb32 <- x^2*T2
  d2lamdb0db1 <- 0
  d2lamdb0db2 <- 0
  d2lamdb0db3 <- 0
  d2lamdb1db2 <- 0
  d2lamdb1db3 <- 0
  d2lamdb2db3 <- 0
  d2lamdb2db3 <- x*T2
  
  
  hess00 <- sum((d2ldp2*dpdb0 + d2ldpdlam*dlamdb0)*dpdb0 + d2pdb02*dldp + (d2ldpdlam*dpdb0 + d2ldlam2*dlamdb0)*dlamdb0 + d2lamdb02*dldlam) 
  hess01 <- sum((d2ldp2*dpdb0 + d2ldpdlam*dlamdb0)*dpdb1 + d2pdb0db1*dldp + (d2ldpdlam*dpdb0 + d2ldlam2*dlamdb0)*dlamdb1 + d2lamdb0db1*dldlam) 
  hess02 <- sum((d2ldp2*dpdb0 + d2ldpdlam*dlamdb0)*dpdb2 + d2pdb0db2*dldp + (d2ldpdlam*dpdb0 + d2ldlam2*dlamdb0)*dlamdb2 + d2lamdb0db2*dldlam) 
  hess03 <- sum((d2ldp2*dpdb0 + d2ldpdlam*dlamdb0)*dpdb3 + d2pdb0db3*dldp + (d2ldpdlam*dpdb0 + d2ldlam2*dlamdb0)*dlamdb3 + d2lamdb0db3*dldlam) 

  hess11 <- sum((d2ldp2*dpdb1 + d2ldpdlam*dlamdb1)*dpdb1 + d2pdb12*dldp + (d2ldpdlam*dpdb1 + d2ldlam2*dlamdb1)*dlamdb1 + d2lamdb12*dldlam) 
  hess12 <- sum((d2ldp2*dpdb1 + d2ldpdlam*dlamdb1)*dpdb2 + d2pdb1db2*dldp + (d2ldpdlam*dpdb1 + d2ldlam2*dlamdb1)*dlamdb2 + d2lamdb1db2*dldlam) 
  hess13 <- sum((d2ldp2*dpdb1 + d2ldpdlam*dlamdb1)*dpdb3 + d2pdb1db3*dldp + (d2ldpdlam*dpdb1 + d2ldlam2*dlamdb1)*dlamdb3 + d2lamdb1db3*dldlam) 

  hess22 <- sum((d2ldp2*dpdb2 + d2ldpdlam*dlamdb2)*dpdb2 + d2pdb22*dldp + (d2ldpdlam*dpdb2 + d2ldlam2*dlamdb2)*dlamdb2 + d2lamdb22*dldlam) 
  hess23 <- sum((d2ldp2*dpdb2 + d2ldpdlam*dlamdb2)*dpdb3 + d2pdb2db3*dldp + (d2ldpdlam*dpdb2 + d2ldlam2*dlamdb2)*dlamdb3 + d2lamdb2db3*dldlam) 

  hess33 <- sum((d2ldp2*dpdb3 + d2ldpdlam*dlamdb3)*dpdb3 + d2pdb32*dldp + (d2ldpdlam*dpdb3 + d2ldlam2*dlamdb3)*dlamdb3 + d2lamdb32*dldlam) 
  
  
  hess <- matrix(c(hess00, hess01, hess02, hess03, 
                   hess01, hess11, hess12, hess13,
                   hess02, hess12, hess22, hess23,
                   hess03, hess13, hess23, hess33), ncol = 4)
  
  return(list(loglik, grad, hess))  
}

#### Data ####
dat<-read.table(file='./GitHub/ExtensiveAssignments/EA 1/data/greenbeandat.txt', header=TRUE)[,c(1,3,4,5)]
dat$date <-  ymd(as.character(dat$date))

#### Sample Stores ####
samp_store<-function(n, seed, dat){
  set.seed(seed)
  s<-sample(unique(dat$store),n)
  dat[dat$store %in% s,]
}

dat1<-samp_store(1, 5, dat=dat)

#### Cheat w/Optim ####

loglik <- function(par, dat){
  b0<-par[1]
  b1<-par[2]
  b2<-par[3]
  b3<-par[4]
  
  x <- dat$price
  y <- dat$mvm
  
  T1 <- exp(b0 + b1*x)
  T2 <- exp(b2 + b3*x)  
  
  p <- T1/(1+T1)
  lambda <- T2
  
  delta<-rep(0,length(y))
  delta[y==0]<-1
  
  sum(delta*log(1-p + p*exp(-lambda)) + (1-delta)*(log(p) - lfactorial(y) + y*log(lambda) - lambda))
}

mles<-optim(c(0,0,0,0), function(u) -loglik(u, dat=dat1), hessian=TRUE)

#### Try ders.newt ####
source("./GitHub/ExtensiveAssignments/EA 1/code/newtraph_kaiser.R")
x0<-c(8, -9, 7, -8)
dat1.newt<-newtraph(ders=ders.bean, dat=dat1, x0=x0)


