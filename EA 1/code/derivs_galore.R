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
  
  delta <- as.numeric(y==0)
  
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
  d2ldlam2 <- delta*p*(1-p)*exp(-lambda)/(1-p + p*exp(-lambda))^2 - (1-delta)*(y/lambda^2)
  d2pdb02 <- (1-T1)*T1/(1+T1)^3
  d2pdb12 <- x^2*(1-T1)*T1/(1+T1)^3
  d2pdb22 <- 0
  d2pdb32 <- 0  
  d2pdb0db1 <- x*(1-T1)*T1/(1+T1)^3
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
dat<-read.table(file='EA 1/data/greenbeandat.txt', header=TRUE)[,c(1,3,4,5)]
dat$date <-  ymd(as.character(dat$date))

#### Sample Stores ####
samp_store<-function(n, seed, dat){
  set.seed(seed)
  s<-sample(unique(dat$store),n)
  dat[dat$store %in% s,]
}

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


#### Try ders.newt ####
source("EA 1/code/newtraph_kaiser.R")
#x0<-c(8, -9, 7, -8)
#dat1.newt<-newtraph(ders=ders.bean, dat=dat1, x0=x0)

library(plyr)
library(reshape2)

est <- dlply(dat, .(store), function(x) {
  x0 <- optim(c(0,0,0,0), function(u) -loglik(u, dat=x))$par
  tryCatch(newtraph(ders=ders.bean, dat=x, x0=x0), 
           error=function(e) list(e))
})

##how many didn't work
sum(ldply(est, function(x) sum(length(x) < 3))[,2])

mles <- ldply(est, function(x) if(length(x) == 3) t(x[[1]]))

ggplot(melt(mles, id.vars="store")) + 
  geom_histogram(aes(x=value, y=..density..)) + 
  facet_wrap(~variable)

ggplot(dat1) + geom_line(aes(x=price, y=p*lambda)) + geom_point(aes(x=price, y=p*lambda)) + 
  geom_line(aes(x=price, y=avg_mvm), data=ddply(dat1, .(price), summarise, avg_mvm=mean(mvm)), colour=I("red")) +
  geom_point(aes(x=price, y=avg_mvm),data=ddply(dat1, .(price), summarise, avg_mvm=mean(mvm)), colour=I("red")) +
  ylab("Marginal Expectation") + xlab("Price") +
  ggtitle(sprintf("Store ID: %d", dat1$store[1]))



plot_by_store <- function(store_id, data=dat) {
  dat1 <- subset(data, store == store_id)
  
  est1 <- est[[as.character(dat1$store[1])]][[1]]
  dat1$p <- exp(est1[1] + est1[2]*dat1$price)/(1 + exp(est1[1] + est1[2]*dat1$price))
  dat1$lambda <- exp(est1[3] + est1[4]*dat1$price)
  
  distns <- ldply(mlply(unique(dat1$price), function(y) mdply(0:max(subset(dat1, price == y)$mvm), function(x) c(mvm=x, price=y, lik=exp(loglik(est1, dat=list(price=y, mvm=x)))))), function(x) x)
  
  gg <- ggplot(dat1) + geom_histogram(aes(x=mvm, y=..density..)) + 
    geom_line(aes(x=mvm, y=lik), data=distns, colour=I("red")) + facet_wrap(~price, scales="free")  +
    ggtitle(sprintf("Store ID: %d", dat1$store[1]))
  
  return(gg)
}

fit_plots <- mlply(1:length(est), function(x) {
  if(length(est[[x]]) == 3) {
    pdf("EA 1/figure/diagnostics_store.pdf")
      plot_by_store(as.numeric(names(est)[x])) 
      dev.off()
  }
})













